from __future__ import absolute_import

import random
from math import log
from operator import itemgetter

import numpy

import Orange

from . import obiGEO
from .obiExpression import ExpressionSignificance_Test

# Utility functiions

log2 = lambda x: log(x, 2.)

def my_ratio(x, y):
    """ compute the log-ratio """
    return log2(x/y)

def sign(x):
    return cmp(x, 0)

def common_domain(data1, data2):
    """Use only attributes that are in both data sets"""
    atts = sorted(set([a.name for a in data1.domain.attributes]).intersection(
               [a.name for a in data2.domain.attributes]))
    new1 = Orange.data.Table(Orange.data.Domain(atts + [data1.domain.classVar], data1.domain), data1)
    new2 = Orange.data.Table(Orange.data.Domain(atts + [data2.domain.classVar], data2.domain), data2)
    return new1, new2

def common_genes(data1, data2, att='gene'):
    common_set = list(set(ex[att].value for ex in data1).intersection(ex[att].value for ex in data2))
    print len(set(ex[att].value for ex in data1))
    print len(common_set)
    return data1.filter(**{att: common_set}), data2.filter(**{att: common_set})

# Normalization

def quantilenorm(data):
    """normalization of microarray data to obtain the same distribution in all chips"""
    chips = data.domain.attributes
    genes = [str(d["gene"]) for d in data]
    d_sort = {}
    for c in chips:
        dc = [(float(d[c]), str(d["gene"])) for d in data]
        dc.sort()
        d_sort[c.name] = dc
    d_sort2 = dict([(str(d["gene"]),{}) for d in data])
    for i in range(len(data)):
        genes_c = [(d_sort[c.name][i][1], c.name) for c in chips]
        mean_row = numpy.mean([d_sort[c.name][i][0] for c in chips])
        for g, d in genes_c:
            d_sort2.get(g, {}).update(dict([(d, mean_row)]))
    data_norm = Orange.data.Table(data.domain)
    for i, d in enumerate(data):
        g = str(d["gene"])
        ex = [d_sort2[g][c.name] for c in chips]
        data_norm.append(ex)
        data_norm[i]["gene"] = g
    return data_norm

def get_mediannorm(modeldata):
    """returns the function for median scale normalization"""
    globalmedian = numpy.median([float(d[c]) for d in modeldata for c in modeldata.domain.attributes])
    def normalize(data):
        chips = data.domain.attributes
        medians = [numpy.median([float(d[c]) for d in data]) for c in chips]
        data_norm = Orange.data.Table(data.domain)
        data_norm.domain.add_metas(data.domain.get_metas())
        for i, d in enumerate(data):
            ex = [(d[c.name].value*globalmedian)/medians[k] for k,c in enumerate(chips)]
            data_norm.append(ex)
            for m in data.domain.get_metas():
                data_norm[i][m] = data[i][m]
        return data_norm
    return normalize

def medianscalenorm(data, newsamples=None, intersection=True):
    """normalization of microarray data to center the distributions in all chips on the same global median"""
    if not newsamples:
        n = get_mediannorm(data)
        return n(data), None
    else:
        if intersection:
            idata, inewsamples = common_genes(data, newsamples)
        n = get_mediannorm(idata)
        return n(idata), n(inewsamples)

def normalize(data1, data2=None, type='median'):
    if type == 'quantile':
        ndata1 = quantilenorm(data1)
        return ndata1, medianscalenorm(ndata1, data2)[1] if data2 else None
    elif type == 'median':
        return medianscalenorm(data1, data2)
    else:
        return Error


# Gene filtering

def compute_diff_pairs(at_a, d):
    """ computes the pairwise differences between the replicates """
    differ = []
    for i in range(len (at_a)-1):
        for j in range(i+1, len(at_a)):
            differ.append(log2(d[at_a[i]])-log2(d[at_a[j]]))
    return differ

def costruct_series(data, attr_set, differences=False):
    """ Constructs the time series by averaging the replicates of the same time point """
    serie = dict([(str(d["gene"]), [numpy.mean([d[at] for at in at_samples])
                    for at_name, at_samples in attr_set]) for d in data])
    if differences:
        """store the differences between replicates while creating the time series"""
        differences = []
        r = [[differences.extend(compute_diff_pairs(at_samples, d))
              for at_name, at_samples in attr_set] for d in data]
        return serie, differences
    else:
        return serie

def costruct_control_series(serie, num_points, control):
    """ Creation of the control series (for each gene) as a constant profile equal to expression at time 0 """
    ##print "control series..."
    serie_0 = {}
    if control == "avg":
        serie_0.update(dict([(g,[numpy.mean(serie[g]) for i in range(num_points)]) for g in serie]))
    elif control == "t0":
        serie_0.update(dict([(g,[serie[g][0] for i in range(num_points)]) for g in serie]))
    return serie_0

def compute_area(vals, t, baseline=0.):
    return sum(a_pair((t[i], vals[i]), (t[i+1], vals[i+1]), baseline) \
        for i in xrange(len(vals)-1))

def a_pair(p1, p2, baseline):
    """Area under the line bounded by a pair of two points (x,y) with
    respect to baseline. Both parts around the diagonal are
    positive."""
    x1,y1 = p1
    x2,y2 = p2
    x2 = x2-x1 #same start
    x1 = 0.0
    a = y1-baseline
    b = y2-baseline
    if a*b >= 0: #both on one side
        return abs(x2*(b+a)/2.0)
    else:
        xp = -a * x2 / float(y2-y1)
        return (abs(xp * a) + abs((x2-xp)*b)) / 2.0

def uniform_time_scale(attr_set):
    """ Obtains time points with a unique measure (the lowest among [min,h,d]) present in data"""
    ord_measures = ["min","h","d"]
    converter = {"min":{"h":60, "d":24*60}, "h": {"d":24}}
    measures = list(set([t.split(" ")[1] for t,s in attr_set]))
    if len(measures) == 1:
        time_points = [float(t.split(" ")[0]) for t,s in attr_set]
    else:
        first_measure = min([(ord_measures.index(m),m) for m in measures])[1]
        time_points = [float(t.split(" ")[0]) for t, s in attr_set
                       if t.split(" ")[1] == first_measure]
        time_points.extend([float(t.split(" ")[0]) * converter[first_measure][t.split(" ")[1]]
                            for t, s in attr_set if t.split(" ")[1] != first_measure])
    return time_points

def AREA(data, attr_set, control='t0', weighted=False, auto=False, perc=99):
    """ AREA (Area Under the Curve) filtering method """
    from matplotlib.mlab import prctile
    if weighted:
        time_points = uniform_time_scale(attr_set)
    else:
        time_points = range(len(attr_set))

    # Monte Carlo approach to create the null distribution of the areas
    if auto: # null distribution
        serie, differences = costruct_series(data, attr_set, auto)
        serie_campionata = {}
        area = []
        for i in range(20000):
            serie_campionata = random.sample(differences, len(time_points))
            area.append(compute_area(serie_campionata, time_points))
        area_threshold = prctile(area, perc)
    else:
        serie = costruct_series(data, attr_set, auto)

    serie_0 = costruct_control_series(serie, len(time_points), control)

    # Gene filtering
    areas = []
    for g in serie:
        diff_s = [log2(serie[g][i]) - log2(serie_0[g][i]) for i in range(len(time_points))]
        area_diff = compute_area(diff_s, time_points);
        if not auto or area_diff > area_threshold:
            areas.append((g, area_diff))
    return areas

class ExpressionSignificance_AREA(ExpressionSignificance_Test):
    def __call__(self, target=None):
        attr_set = {}
        for a in self.data.domain.attributes:
            attr_set[a.attributes['time']] = attr_set.get(a.attributes['time'], []) + [a.name]
        scores = AREA(self.data, sorted(attr_set.items()))
        gene2ind = dict((g, i) for i,g in enumerate(ex['gene'].value for ex in self.data))
        return [(gene2ind[g], s) for g, s in scores]

def FC(data, attr_set, control='t0', thr=2, auto=False, p_thr=0.2):
    """ Gene filtering based on the number of FC of all time points with the control series > thr """
    serie = costruct_series(data, attr_set, False)
    num_points = len(serie.values()[0])
    serie_0 = costruct_control_series(serie, num_points, control)
    fc = [(g, len([0 for i in range(num_points) if abs(my_ratio(s[i], serie_0[g][i])) >= thr]))
          for g, s in serie.items()]
    if auto:
        thr_points = round(p_thr*num_points)
        fc = [(g, v) for g, v in fc if v >= thr_points]
    return fc

class ExpressionSignificance_FCts(ExpressionSignificance_Test):
    def __call__(self, target=None):
        attr_set = {}
        for a in self.data.domain.attributes:
            attr_set[a.attributes['time']] = attr_set.get(a.attributes['time'], []) + [a.name]
        scores = FC(self.data, sorted(attr_set.items()))
        return zip(self.keys, map(itemgetter(1), scores))

def spearmanr_filter(data, limit=1000):
    """ Spearman ranks gene filtering """
    from scipy.stats import spearmanr
    time = [a.attributes['time'] for a in data.domain.attributes]
    exs = sorted(data, reverse=True, key=lambda ex: abs(spearmanr([a.value for a in ex], time)[0]))
    return [str(ex['gene']) for ex in exs[:limit]]


def signed_PCA(data):
    pca = Orange.projection.linear.PCA(data, standardize=False, max_components=1)
    classifier = lambda X: [x[0].value for x in pca(X)]
    predictions = classifier(data)
    classes = [ex.getclass().value for ex in data]
    n = 0
    for i1,c1 in enumerate(classes):
        for i2,c2 in enumerate(classes[:i1]):
            n += cmp(c1,c2) * cmp(predictions[i1], predictions[i2])
    if n < 0:
        def invert(X):
            y = classifier(X)
            return -y if type(y) == float else [-x for x in y]
        return invert
    else:
        return classifier

signed_PCA.name = 'PCA'


def conttime(data, d):
    for a in data.domain.attributes:
        a.attributes['time'] = d[a.attributes['time']]

def conv(attr_set, ticks=True):
    """Obtain time points with a unique measure (the lowest among [min,h,d]) present in data"""
    ord_measures = ["min","h","d"]
    converter = {"min":{"h":60, "d":24*60}, "h": {"d":24}}
    measures = list(set([t.split(" ")[1] for t in attr_set]))
    if len(measures) == 1:
        time_points = [(t, float(t.split(" ")[0])) for t in attr_set]
    else:
        first_measure = min([(ord_measures.index(m),m) for m in measures])[1]
        time_points = [(t, float(t.split(" ")[0])) for t in attr_set if t.split(" ")[1] == first_measure]
        time_points.extend([(t, float(t.split(" ")[0]) * converter[first_measure][t.split(" ")[1]])
                            for t in attr_set if t.split(" ")[1] != first_measure])
    time_points.sort(key=itemgetter(1))
    if ticks:
        time_points = [(t[0],float(i)) for i,t in enumerate(time_points)]
    return dict(time_points)


def get_projections(data1, data2=None):
    labels1 = list(a.attributes['time'] for a in data1.domain.attributes)
    tdata1 = obiGEO.transpose(data1)
    if data2:
        labels2 = list('[%s]' % a.attributes['time'] for a in data2.domain.attributes)
        tdata2 = obiGEO.transpose(data2)
        tdata1, tdata2 = common_domain(tdata1, tdata2)
        classifier = signed_PCA(tdata1)
        proj1 = classifier(tdata1)
        proj2 = classifier(tdata2)
    else:
        classifier = signed_PCA(tdata1)
        proj1 = classifier(tdata1)
        proj2, labels2 = [], []
    return proj1, labels1, proj2, labels2


############

if False and __name__ == '__main__':
    from . import obiGEO
    # Data set 1
    data1 = obiGEO.GDS('GDS2666').getdata(report_genes=True, transpose=False)
    labels1 = list(a.attributes['time'] for a in data1.domain.attributes)
    attr_set = list(set(a.attributes['time'] for a in data1.domain.attributes))
    convd = conv(attr_set)
    conttime(data1, convd)
    # Data set 2
    data2 = obiGEO.GDS('GDS2667').getdata(report_genes=True, transpose=False)
    labels2 = list(a.attributes['time'] for a in data2.domain.attributes)
    attr_set = list(set(a.attributes['time'] for a in data2.domain.attributes))
    convd = conv(attr_set)
    conttime(data2, convd)
    # Normalize data set 1
    ndata1, _ = normalize(data1, type='quantile')
    # Filtering
    attr_set = {}
    for a in ndata1.domain.attributes:
        attr_set[a.attributes['time']] = attr_set.get(a.attributes['time'], []) + [a.name]
    scores = AREA(ndata1, sorted(attr_set.items()))
    genes = map(itemgetter(0), sorted(scores, key=itemgetter(1))[:1000])
    fdata1 = ndata1.filter(gene=genes)
    # Rescale data set 2
    ttrain, ttest = normalize(fdata1, data2)
    # Model construction and prediction
    # train = obiGEO.transpose(ttrain)
    # test = obiGEO.transpose(ttest)
    # cdtrain, cdtest = common_domain(train, test)
    # classifier = signed_PCA(cdtrain)
    # proj1 = classifier(cdtrain)
    # proj2 = classifier(cdtest)



