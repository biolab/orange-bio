# compatibility with Orange 2 and 3

import Orange

OR3 = True if Orange.__version__ >= "3" else False

if OR3:
    from Orange.data import DiscreteVariable, ContinuousVariable, StringVariable
else:
    from Orange.feature import Discrete as DiscreteVariable
    from Orange.feature import Continuous as ContinuousVariable
    from Orange.feature import String as StringVariable

if OR3:
    import numpy
    unknown = NAN = float("nan")
    isunknown = numpy.isnan
else:
    unknown = NAN = '?'
    isunknown = lambda v: v == NAN


def create_domain(at, cl, metas):
    if OR3:
        return Orange.data.Domain(at, cl, metas=metas)
    else:
        domain  = Orange.data.Domain(at, cl)
        if metas:
            if isinstance(metas, dict):
                metas = sorted(metas.items())
            else:
                metas = zip([ StringVariable.new_meta_id() for _ in metas ], metas)
            domain.add_metas(dict((StringVariable.new_meta_id(), ma) for mi, ma in metas))
        return domain


def create_table(domain, X, Y, metas):
    classvar = domain.class_var
    if OR3:
        metaatts = domain.metas
        if Y:
            Y = numpy.array([[classvar.to_val(row)] for row in Y],
                            dtype=float)
        if metas:
            metas = numpy.array([[c.to_val(v) for c, v in zip(metaatts, row)]
                                 for row in metas],
                                dtype=object)
        data = Orange.data.Table(domain, numpy.asarray(X), Y=Y, metas=metas)
    else:
        metaatts = [ at for i,at in sorted(domain.getmetas().items()) ]
        insts = []
        for i in range(len(X)):
            if Y:
                inst = Orange.data.Instance(domain, X[i] + [ Y[i] ])
            else:
                inst = Orange.data.Instance(domain, X[i])
            if metaatts:
                for ma,mv in zip(metaatts, metas[i]):
                    inst[ma] = mv
            insts.append(inst)
        data = Orange.data.Table(domain, insts)
    return data

def get_metas(domain):
    if OR3:
        return domain.metas
    else:
        return [at for i, at in sorted(domain.getmetas().items())]
