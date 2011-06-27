import random
from math import log

import numpy
from matplotlib.mlab import prctile

import Orange
import obiGEO


# Normalization

def quantilenorm(data):
    """ normalization of microarray data to obtain the same distribution in all chips """
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
    """ returns the function for median scale normalization """
    globalmedian = numpy.median([float(d[c]) for d in modeldata for c in modeldata.domain.attributes])
    def normalize(data):
        chips = data.domain.attributes
        medians = [numpy.median([float(d[c]) for d in data]) for c in chips]
        data_norm = Orange.data.Table(data.domain)
        for i, d in enumerate(data):
            g = str(d["gene"])
            ex = [(d[c.name].value*globalmedian)/medians[k] for k,c in enumerate(chips)]
            data_norm.append(ex)
            data_norm[i]["gene"] = g
        return data_norm
    return normalize

def medianscalenorm(data, newsamples=None):
    """ normalization of microarray data to obtain the distributions in all chips centered on the same global median """
    n = get_mediannorm(data)
    if not newsamples:
        return n(data), None
    else:
        return n(data), n(newsamples)

def normalize(data1, data2=None, type='median'):
    if type == 'quantile':
        ndata1 = quantilenorm(data1)
        return ndata1, medianscalenorm(ndata1, data2)[1] if data2 else None
    elif type == 'median':
        return medianscalenorm(data1, data2)
    else:
        return Error


# Gene filtering

log2 = lambda x: log(x, 2.)

def my_ratio(x, y):
    """ compute the log-ratio """
    return log2(x/y)

def sign(int):
    """ returns the sign of int """
    if(int < 0):
        return -1;
    elif(int > 0):
        return 1;
    else:
        return 0;

def compute_diff_pairs(at_a, d):
    """ computes the pairwise differences between the replicates """
    differ = []
    for i in range(len (at_a)-1):
        for j in range(i+1, len(at_a)):
            differ.append(log2(d[at_a[i]])-log2(d[at_a[j]]))
    return differ

def costruct_series(data, attr_set, f):
    """ Constructs the time series by averaging the replicates of the same time point """
    ##print "creating the time series..."
    serie = dict([(str(d["gene"]), [numpy.mean([d[at] for at in at_samples]) for at_name, at_samples in attr_set]) for d in data])
    if f == "AREA":
        """ for AREA filtering, store the differences between replicates while creating the time series """
        differences = []
        r = [[differences.extend(compute_diff_pairs(at_samples, d)) for at_name, at_samples in attr_set] for d in data]
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

def compute_area(serie, tempi):
    """ Function that calculates the area between 2 time series """
    # if the ith point and the next point have the same sign it computes the area by the
    # integral with the trapezoidal method, else it computes the area of the two triangles found with linear interpolation
    area_pieces = []
    for i in range(len(serie)-1):
        if sign(serie[i]) == sign(serie[i+1]):
            area_pieces.append(abs(numpy.trapz([serie[i], serie[i+1]], [tempi[i], tempi[i+1]])))
        else:
            area_triangle_1 = float(serie[i]*(tempi[i+1]-tempi[i]))/(serie[i]-serie[i+1])*abs(serie[i])/2
            area_triangle_2 = float(serie[i+1]*(tempi[i]-tempi[i+1]))/(serie[i]-serie[i+1])*abs(serie[i+1])/2
            area_pieces.append(area_triangle_1 + area_triangle_2)

    area = sum(area_pieces);
    return area

def uniform_time_scale(attr_set):
    """ Obtains time points with a unique measure (the lowest among [min,h,d]) present in data"""
    ord_measures = ["min","h","d"]
    converter = {"min":{"h":60, "d":24*60}, "h": {"d":24}}
    measures = list(set([t.split(" ")[1] for t,s in attr_set]))
    if len(measures) == 1:
        time_points = [float(t.split(" ")[0]) for t,s in attr_set]
    else:
        first_measure = min([(ord_measures.index(m),m) for m in measures])[1]
        time_points = [float(t.split(" ")[0]) for t,s in attr_set if t.split(" ")[1] == first_measure]
        time_points.extend([float(t.split(" ")[0])*converter[first_measure][t.split(" ")[1]] for t,s in attr_set if t.split(" ")[1] != first_measure])
    return time_points

def AREA(data, attr_set, perc, control, scale, limit):
    """ AREA (Area Under the Curve) filtering method """
    if scale == "weighted":
        time_points = uniform_time_scale(attr_set)
    else:
        time_points = range(len(attr_set))

    serie, differences = costruct_series(data, attr_set, "AREA")
    serie_0 = costruct_control_series(serie, len(time_points), control)

    """ Monte Carlo approach to create the null distribution of the areas """
    ##print "null distribution..."
    serie_campionata = {}
    area = []
    for i in range(20000):
        serie_campionata = random.sample(differences, len(time_points));
        area.append(compute_area(serie_campionata, time_points));

    area_soglia = prctile(area,perc);


    """ Gene filtering based on the threshold (selected percentile of null distribution of areas) """
    ##print "filtering..."
    diff_expr = []
    areas = []
    for g in serie:
        diff_s = [log2(serie[g][i]) - log2(serie_0[g][i]) for i in range(len(time_points))]
        area_diff = compute_area(diff_s, time_points);
        areas.append((area_diff,g))
        if area_diff > area_soglia:
            diff_expr.append(g)
    areas.sort(reverse = True)
    if limit == 0:
        de_genes = diff_expr
    else:
        de_genes = [d[1] for d in areas][:limit]

    ##print "%d genes differentially expressed" %len(de_genes)
    return de_genes

def FC(data, attr_set, thr, p_thr, control, limit):
    """ Gene filtering based on the number of FC of all time points with the control series > thr """
    serie = costruct_series(data, attr_set, "fc")
    num_points = len(serie.values()[0])
    serie_0 = costruct_control_series(serie, num_points, control)
    thr_points = round(p_thr*num_points)
    ##print "filtering..."
    if limit == 0:
        de_genes = [gene for gene,s in serie.items() if len([True for i in range(num_points) if abs(my_ratio(s[i], serie_0[gene][i])) >= log2(thr)]) >= thr_points]
    else:
        fc = [(len([True for i in range(num_points) if abs(my_ratio(s[i], serie_0[gene][i])) >= log2(thr)]),gene) for gene,s in serie.items()]
        fc.sort(reverse = True)
        de_genes = [d[1] for d in fc][:limit]

    ##print "%d genes differentially expressed" %len(de_genes)
    return de_genes

def spearmanr_filter(data, limit=1000):
    """ Spearman ranks gene filtering """
    from scipy.stats import spearmanr
    time = [a.attributes['time'] for a in data.domain.attributes]
    exs = sorted(data, reverse=True, key=lambda ex: abs(spearmanr([a.value for a in ex], time)[0]))
    return [str(ex['gene']) for ex in exs[:limit]]


