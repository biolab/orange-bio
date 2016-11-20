"""GenAPI"""
import json
import io
import gzip
import requests


from genesis import Genesis, GenData
from Orange.data import ContinuousVariable, StringVariable, Domain, Table, DiscreteVariable
from .tools import transpose_table
from .. import dicty

CallBack = dicty.CallBack
transformValues = dicty.transformValues
averageAttributes = dicty.averageAttributes
example_tables = dicty.example_tables

view_model = ['experiment', 'growth', 'genotype', 'treatment', 'strain', 'time', 'replicate']
DEFAULT_EMAIL = 'anonymous@genialis.com'
DEFAULT_PASSWD = 'anonymous'
DEFAULT_URL = 'https://dictyexpress.research.bcm.edu'


class GenAPI(object):
    """
    Python module that leverages Genesis PyAPI (Python API for accsess to DictyExpress database).
    It supports connection to the server and data retrieval functionalities.
    """
    def __init__(self, email=DEFAULT_EMAIL, password=DEFAULT_PASSWD, url=DEFAULT_URL):

        self._gen = Genesis(email, password, url)
        self.email = email

    def fetch_etc_objects(self, callback=lambda: None):
        """ Function downloads all available :obj:`GenData` etc objects from DictyExpress database.

        Returns:
            :obj:`list`: of :obj:`GenData` objects

        """
        cbc = CallBack(1, callback)
        try:
            list_of_experiments = self._gen.api.data.get(case_ids__contains='5535115cfad58d5e03006217', status='done',
                                                         type__startswith='data:etc:')['objects']

        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        cbc()
        store_experiments = [GenData(exp, self._gen) for exp in list_of_experiments]
        cbc.end()
        return store_experiments

    def download_etc_data(self, gen_data_id, callback=lambda: None):
        """ Function downloads etc data of a chosen experiment from the server.

        Args:
            gen_data_id (str): id of :obj:`GenData` object to download.

        Returns:
             :obj:`dict`: data in json like format


        """
        cbc = CallBack(3, callback, callbacks=70)

        try:
            response = next(self._gen.download([gen_data_id], 'output.etcfile'))
        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        cbc()
        if not response.ok:
            response.raise_for_status()

        response_gzipped = io.BytesIO(response.content)
        cbc()
        response_content = io.TextIOWrapper(gzip.GzipFile(fileobj=response_gzipped), encoding="utf-8")
        cbc()

        try:
            json_object = json.load(response_content)
        except ValueError as e:
            raise ValueError('Downloaded data is not a valid JSON') from e

        cbc.end()
        return json_object

    def etc_to_table(self, etc_json, time_var=False, callback=lambda: None):
        """ Converts data from Json to :obj:`Orange.data.table`

        Args:
            etc_json (dict): Data in json like format
            time_var (bool): Create column of time points. Default is set to False.
        Returns:
            :obj:`Orange.data.Table`
        """
        cbc = CallBack(2, callback, callbacks=30)

        variables = []
        time_point = 1
        for time in etc_json['etc']['timePoints']:
            var = ContinuousVariable('TP ' + str(time_point))
            var.attributes['Time'] = str(time)
            variables.append(var)
            time_point += 1

        meta_attr = StringVariable.make('Gene')
        domain = Domain(variables, metas=[meta_attr])
        cbc()

        table = []
        for row in etc_json['etc']['genes']:
            gene_expression = [exp for exp in etc_json['etc']['genes'][row]]
            gene_expression.append(row)
            table.append(gene_expression)

        orange_table = Table(domain, table)

        if time_var:
            orange_table = transpose_table(orange_table)
            cbc()

        cbc.end()
        return orange_table
