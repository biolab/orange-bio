"""Tools"""
from .. import dicty
from Orange.data import ContinuousVariable, StringVariable, TimeVariable, Domain, Table


CallBack = dicty.CallBack
transformValues = dicty.transformValues
averageAttributes = dicty.averageAttributes
example_tables = dicty.example_tables


def transpose_table(table):
    """
    Transpose the rows and columns of the table.

    Args:
        table: Data in :obj:`Orange.data.Table`

    Returns:
         Transposed :obj:`Orange.data.Table`. (Genes as columns)
    """
    attrs = table.domain.attributes
    attr = [ContinuousVariable.make(ex['Gene'].value) for ex in table]
    #  Set metas
    new_metas = [StringVariable.make(name) if name is not 'Time' else TimeVariable.make(name)
                 for name in sorted(table.domain.variables[0].attributes.keys())]
    domain = Domain(attr, metas=new_metas)
    meta_values = [[exp.attributes[var.name] for var in domain.metas] for exp in attrs]

    return Table(domain, table.X.transpose(), metas=meta_values)
