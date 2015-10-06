import os
import sqlite3
import tarfile
import shutil
import tempfile
import collections
import textwrap

from collections import namedtuple

try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import six

__all__ = ["Taxonomy"]

pjoin = os.path.join

TAXDUMP_URL = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

_taxon = namedtuple(
    "taxon",
    ["tax_id", "parent_tax_id", "name", "synonyms", "rank"]
)


def namedtuple_repr_pretty(self, p, cycle):
    name = type(self).__name__
    if cycle:
        p.text("{0}(...)".format("name"))
    else:
        with p.group(len(name) + 1, "{0}(".format(name), ")"):
            for field, val in zip(self._fields, self):
                p.text(field + "=")
                p.pretty(val)
                if field != self._fields[-1]:
                    p.text(",")
                    p.breakable()


class taxon(_taxon):
    _repr_pretty_ = namedtuple_repr_pretty


class Taxonomy(collections.Mapping):
    SCHEMA_VERSION = (0, 0, 1)

    def __init__(self, taxdb):
        self._con = sqlite3.connect(taxdb, timeout=15)
        self._con.execute("""
            CREATE INDEX IF NOT EXISTS
                index_names_tax_id ON names(tax_id)
        """)

    def __node_query(self, tax_id):
        c = self._con.execute("""
            SELECT nodes.tax_id, nodes.parent_tax_id, ranks.rank
            FROM nodes INNER JOIN ranks USING(rank_id)
            WHERE tax_id = ?
        """, (tax_id,))
        if not c:
            raise KeyError(tax_id)
        else:
            return next(c)

    def __getitem__(self, tax_id):
        if not isinstance(tax_id, six.string_types):
            raise TypeError("expected a string, got {}"
                            .format(type(tax_id).__name__))

        c = self._con.execute("""
            SELECT names.name, name_classes.name_class
            FROM names INNER JOIN name_classes USING (name_class_id)
            WHERE names.tax_id = ?
        """, (tax_id,))

        names = list(c)
        if not names:
            raise KeyError(tax_id)

        scientific_name = None
        for name, name_class in names:
            if name_class == "scientific name":
                scientific_name = name
        node = self.__node_query(tax_id)
        _, parent, rank = node
        return taxon(str(tax_id), str(parent), scientific_name, names, rank)

    def __iter__(self):
        c = self._con.execute("SELECT tax_id FROM nodes")
        return (str(r[0]) for r in c)

    def __len__(self):
        c = self._con.execute("SELECT COUNT(*) FROM nodes")
        return next(c)[0]

    def search(self, name, exact=True):
        # First ensure the name column is indexed.
        self._con.execute("""
            CREATE INDEX IF NOT EXISTS
                index_names_name ON names(name)""")

        if not exact:
            name = "{0}%".format(name)
            operator = "LIKE"
        else:
            operator = "="

        c = self._con.execute("""
            SELECT DISTINCT(tax_id)
            FROM names
            WHERE names.name {operator} ?
            """.format(operator=operator), (name,))
        return (str(r[0]) for r in c)

    def lineage(self, tax_id):
        lineage = []
        while True:
            parent = self.parent_tax_id(tax_id)
            if parent is not None:
                lineage.append(parent)
                tax_id = parent
            else:
                break
        return list(reversed(lineage))

    def parent_tax_id(self, tax_id):
        if not isinstance(tax_id, six.string_types):
            raise TypeError("Expected a string")

        node = self.__node_query(tax_id)
        _, parent, _ = node
        if parent == int(tax_id):
            # The root (with tax_id 1) has itself as a parent
            return None
        else:
            return str(parent)

    def child_tax_ids(self, tax_id):
        if not isinstance(tax_id, six.string_types):
            raise TypeError("Expected a string")

        c = self._con.execute("""
            SELECT tax_id
            FROM nodes
            WHERE parent_tax_id = ?
        """, (tax_id,))
        children = list(str(r[0]) for r in c if r[0] != 1)
        return children

    def name(self, tax_id):
        return self[tax_id].name

    def synonyms(self, tax_id):
        return self[tax_id].synonyms

    @classmethod
    def initialize(cls, db_filename, taxdump=None):

        tempd = None
        if taxdump is None:
            tempd = tempfile.mkdtemp()
            cls.download(tempd)
            taxdump = pjoin(tempd, "taxdump.tar.gz")

        try:
            cls.init_db(db_filename, tarfile.open(taxdump))
        finally:
            if tempd is not None:
                shutil.rmtree(tempd)

    @classmethod
    def download(cls, download_dir):
        """
        Download the taxonomy archive from the ncbi ftp server.

        :param str download_dir: Download directory.

        """
        stream = urlopen(TAXDUMP_URL, timeout=30)

        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)

        with open(pjoin(download_dir, "taxdump.tar.gz"), "wb") as f:
            shutil.copyfileobj(stream, f)

    @classmethod
    def init_db(cls, dbfilename, taxdump):

        con = sqlite3.connect(dbfilename)
        cursor = con.cursor()

        indexes = ["index_names_tax_id",
                   "index_names_name"]

        for index in indexes:
            cursor.execute("DROP INDEX IF EXISTS %s" % index)

        for table in ["nodes", "name_classes", "names", "ranks"]:
            cursor.execute("DROP TABLE IF EXISTS %s" % table)

        script = textwrap.dedent("""
            CREATE TABLE ranks (
                rank_id INTEGER PRIMARY KEY ASC,
                rank TEXT UNIQUE
            );

            CREATE TABLE name_classes (
                name_class_id INTEGER PRIMARY KEY ASC,
                name_class TEXT UNIQUE
            );

            CREATE TABLE nodes (
                tax_id INTEGER PRIMARY KEY ASC,
                parent_tax_id INTEGER,
                rank_id INTEGER REFERENCES ranks(rank_id)
            );

            CREATE TABLE names (
                tax_id INTEGER REFERENCES nodes(tax_id),
                name text COLLATE NOCASE,
                name_class_id INTEGER REFERENCES name_classes(name_class_id)
            );
        """)

        cursor.executescript(script)

        nodes = taxdump.extractfile("nodes.dmp")

        def iter_rows(lines):
            for line in lines:
                if not line.strip():
                    continue
                yield tuple(line.decode("utf-8").rstrip("\t\n|").split("\t|\t"))

        # take the tax_id, parent_tax_id, rank entries
        nodes = [row[:3] for row in iter_rows(nodes)]

        ranks = list(enumerate(sorted(set(rank for _, _, rank in nodes))))
        rank_id = {rank: i for i, rank in ranks}

        cursor.executemany("INSERT INTO ranks VALUES (?, ?)", ranks)

        cursor.executemany("INSERT INTO nodes VALUES (?, ?, ?)",
                           ((int(tax_id), int(parent_tax_id), rank_id[rank])
                            for tax_id, parent_tax_id, rank in nodes))

        names = taxdump.extractfile("names.dmp")

        def iter_names(rows):
            for tax_id, name, unique_name, name_class in rows:
                if unique_name:
                    name = unique_name
                yield tax_id, name, name_class

        names = list(iter_names(iter_rows(names)))
        name_classes = list(enumerate(set(name_class
                                          for _, _, name_class in names)))
        name_class_id = {name_class: i for i, name_class in name_classes}

        cursor.executemany("INSERT INTO name_classes VALUES (?, ?)",
                           name_classes)
        cursor.executemany("INSERT INTO names VALUES (?, ?, ?)",
                           ((int(tax_id), name, name_class_id[name_class])
                            for tax_id, name, name_class in names))

        con.commit()
        con.close()
