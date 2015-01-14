from __future__ import absolute_import

import os
import urllib2
import shutil
import gzip
import zipfile
import xml.dom.minidom as minidom
import csv
import sqlite3
import errno
import posixpath
import textwrap

from StringIO import StringIO
from collections import defaultdict, namedtuple
from operator import itemgetter

from Orange.utils import serverfiles
from Orange.utils import ConsoleProgressBar, wget

from . import taxonomy

from .taxonomy import pickled_cache


def mkdir_p(path, mode=0o777):
    try:
        os.makedirs(path, mode)
    except OSError as err:
        if err.errno == errno.EEXIST:
            pass
        else:
            raise


class PPIDatabase(object):
    """
    A general interface for protein-protein interaction database access.

    An example::

        >>> ppidb = MySuperPPIDatabase()
        >>> ppidb.organisms() # List all organisms (taxids)
        ['...

        >>> ppidb.ids() # List all protein ids
        ['...

        >>> ppidb.ids(taxid="9606") # List all human protein ids.
        ['...

        >>> ppidb.links() # List all links
        [('...

    """
    def __init__(self):
        pass

    def organisms(self):
        """
        Return all organism ncbi taxonomy ids contained in this database.
        """
        raise NotImplementedError

    def ids(self, taxid=None):
        """
        Return a list of all protein ids. If `taxid` (as returned by
        `organisms()`) is not ``None`` limit the results to ids to
        this organism only.

        """
        raise NotImplementedError

    def synonyms(self, id):
        """
        Return a list of synonyms for primary `id` (as returned by `ids`).
        """
        raise NotImplementedError

    def all_edges(self, taxid=None):
        """
        Return a list of all edges. If `taxid` is not ``None`` return the
        edges for this organism only.

        """
        raise NotImplementedError

    def edges(self, id1, id2=None):
        """
        Return a list of all edges (a list of 3-tuples (id1, id2, score)).
        """
        raise NotImplementedError

    def all_edges_annotated(self, taxid=None):
        """
        Return a list of all edges annotated. If taxid is not None
        return the edges for this organism only.

        """
        res = []
        for id in self.ids(taxid):
            res.extend(self.edges_annotated(id))
        return res

    def edges_annotated(self, id=None):
        """
        Return a list of all edges annotated.
        """
        raise NotImplementedError

    def search_id(self, name, taxid=None):
        """
        Search the database for protein name. Return a list of matching
        primary ids. Use `taxid` to limit the results to a single organism.

        """
        raise NotImplementedError

    def extract_network(self, ids):
        """
        """
        from Orange import network

        graph = network.Graph()
        for id in ids:
            graph.add_node(id, synonyms=",".join(self.synonyms(id)))

        for id in ids:
            edges = self.edges(id)

            for id1, id2, score in edges:
                graph.add_edge(id1, id2, weight=score)

        return graph

    @classmethod
    def download_data(self):
        """
        Download the latest PPI data for local work.
        """
        raise NotImplementedError


class BioGRID(PPIDatabase):
    """
    Access `BioGRID <http://thebiogrid.org>`_ PPI data.

    Example ::

        >>> biogrid = BioGRID()
        >>> print biogrid.organism() # Print a list of all organism ncbi taxis in BioGRID
        [u'10090',...

        >>> print biogrid.ids(taxid="9606") # Print a set of all human protein ids
        [u'110004'

        >>> print biogrid.synonyms("110004") # Print a list of all synonyms for protein id '110004' as reported by BioGRID
        [u'3803', u'CU464060.2', u'CD158b', u'p58.2', u'CD158B1', u'NKAT6']

        >>>

    """

    SCHEMA = [
        ("links",
         """\
         biogrid_interaction_id text,
         biogrid_id_interactor_a text,
         biogrid_id_interactor_b text,
         experimental_system text,
         experimental_system_type text,
         author text,
         pubmed_id text,
         throughput text,
         score real,
         modification text,
         phenotypes text,
         qualifications text,
         tags text,
         source_database text
         """),
        ("proteins",
         """\
         biogrid_id_interactor text,
         entrez_gene_interactor text,
         systematic_name_interactor text,
         official_symbol_interactor text,
         synonyms_interactor text,
         organism_interactor text,
         """)
    ]

    VERSION = "2.0"

    # All column names in the BioGRID tab2 source table.
    FIELDS = ['biogrid_interaction_id',
              'entrez_gene_interactor_a',
              'entrez_gene_interactor_b',
              'biogrid_id_interactor_a',
              'biogrid_id_interactor_b',
              'systematic_name_interactor_a',
              'systematic_name_interactor_b',
              'official_symbol_interactor_a',
              'official_symbol_interactor_b',
              'synonyms_interactor_a',
              'synonyms_interactor_b',
              'experimental_system',
              'experimental_system_type',
              'author',
              'pubmed_id',
              'organism_interactor_a',
              'organism_interactor_b',
              'throughput',
              'score',
              'modification',
              'phenotypes',
              'qualifications',
              'tags',
              'source_database'
              ]

    DOMAIN = "PPI"
    SERVER_FILE = "BIOGRID-ALL.sqlite"

    TAXID_MAP = {
        "4530": "39947",
        "4932": "559292",
        "562": "511145",
        "5833": "36329",
        "2104": None,
        "4754": None,
        "31033": None
    }

    def __init__(self):
        self.filename = serverfiles.localpath_download(
            self.DOMAIN, self.SERVER_FILE)

        # assert version matches
        self.db = sqlite3.connect(self.filename)
        self.init_db_index()

    def organisms(self):
        cur = self.db.execute("select distinct organism_interactor \n"
                              "from proteins")
        return map(itemgetter(0), cur.fetchall())

    def ids(self, taxid=None):
        """
        Return a list of all protein ids (biogrid_id_interactors).
        If `taxid` is not None limit the results to ids from this organism
        only.

        """
        if taxid is None:
            cur = self.db.execute("""\
                select biogrid_id_interactor
                from proteins""")
        else:
            cur = self.db.execute("""\
                select biogrid_id_interactor
                from proteins
                where organism_interactor=?""",
                (taxid,))

        return [t[0] for t in cur.fetchall()]

    def synonyms(self, id):
        """
        Return a list of synonyms for primary `id`.
        """
        cur = self.db.execute("""\
            select entrez_gene_interactor,
                   systematic_name_interactor,
                   official_symbol_interactor,
                   synonyms_interactor
            from proteins
            where biogrid_id_interactor=?""",
            (id,))
        rec = cur.fetchone()
        if rec:
            synonyms = list(rec[:-1]) + \
                       (rec[-1].split("|") if rec[-1] is not None else [])
            return [s for s in synonyms if s is not None]
        else:
            return []

    def all_edges(self, taxid=None):
        """
        Return a list of all edges. If taxid is not None return the
        edges for this organism only.

        """
        if taxid is not None:
            cur = self.db.execute("""\
                select biogrid_id_interactor_a, biogrid_id_interactor_a, score
                from links left join proteins on
                    biogrid_id_interactor_a=biogrid_id_interactor or
                    biogrid_id_interactor_b=biogrid_id_interactor
                where organism_interactor=?
            """, (taxid,))
        else:
            cur = self.db.execute("""\
                select biogrid_id_interactor_a, biogrid_id_interactor_a, score
                from links
            """)
        edges = cur.fetchall()
        return edges

    def edges(self, id):
        """
        Return a list of all interactions where id is a participant
        (a list of 3-tuples (id_a, id_b, score)).

        """

        cur = self.db.execute("""\
            select biogrid_id_interactor_a, biogrid_id_interactor_b, score
            from links
            where biogrid_id_interactor_a=? or biogrid_id_interactor_b=?
        """, (id, id))
        return cur.fetchall()

    def all_edges_annotated(self, taxid=None):
        """
        Return a list of all edges annotated. If taxid is not None
        return the edges for this organism only.

        """
        if taxid is not None:
            cur = self.db.execute("""\
                select *
                from links left join proteins on
                    biogrid_id_interactor_a=biogrid_id_interactor or
                    biogrid_id_interactor_b=biogrid_id_interactor
                where organism_interactor=?
            """, (taxid,))
        else:
            cur = self.db.execute("""\
                select *
                from links
            """)
        edges = cur.fetchall()
        return edges

    def edges_annotated(self, id):
        """ Return a list of all links
        """
        cur = self.db.execute("""\
            select *
            from links
            where biogrid_id_interactor_a=? or biogrid_id_interactor_b=?
        """, (id, id))
        return cur.fetchall()

    def search_id(self, name, taxid=None):
        """
        Search the database for protein name. Return a list of matching
        primary ids. Use `taxid` to limit the results to a single organism.

        """
        # TODO: synonyms_interactor can contain multiple synonyms
        # (should create an indexed table of synonyms)
        if taxid is None:
            cur = self.db.execute("""\
                select biogrid_id_interactor
                from proteins
                where (biogrid_id_interactor=? or
                       entrez_gene_interactor=? or
                       systematic_name_interactor=? or
                       official_symbol_interactor=? or
                       synonyms_interactor=?)
            """, ((name,) * 5))
        else:
            cur = self.db.execute("""\
                select biogrid_id_interactor
                from proteins
                where (biogrid_id_interactor=? or
                       entrez_gene_interactor=? or
                       systematic_name_interactor=? or
                       official_symbol_interactor=? or
                       synonyms_interactor=?)
                      and organism_interactor=?
            """, ((name,) * 5) + (taxid,))
        res = map(itemgetter(0), cur)
        return res

    @classmethod
    def download_data(cls, address):
        """
        Pass the address of the latest BIOGRID-ALL release (in tab2 format).
        """
        stream = urllib2.urlopen(address)
        stream = StringIO(stream.read())
        zfile = zipfile.ZipFile(stream)
        # Expecting only one file.
        filename = zfile.namelist()[0]

        filepath = serverfiles.localpath("PPI", "BIOGRID-ALL.tab2")
        mkdir_p(os.path.dirname(filepath))

        with open(filepath, "wb") as f:
            shutil.copyfileobj(zfile.open(filename, "r"), f)

        cls.init_db(filepath)

    @classmethod
    def init_db(cls, filepath):
        """
        Initialize the sqlite data base from a BIOGRID-ALL.*tab2.txt file
        format.

        """
        dirname = os.path.dirname(filepath)
        rows = csv.reader(open(filepath, "rb"), delimiter="\t")
        rows.next()  # read the header line

        con = sqlite3.connect(os.path.join(dirname, BioGRID.SERVER_FILE))
        con.execute("drop table if exists links")  # Drop old table
        con.execute("drop table if exists proteins")  # Drop old table

        con.execute("""\
            create table links (
                biogrid_interaction_id text,
                biogrid_id_interactor_a text,
                biogrid_id_interactor_b text,
                experimental_system text,
                experimental_system_type text,
                author text,
                pubmed_id text,
                throughput text,
                score real,
                modification text,
                phenotypes text,
                qualifications text,
                tags text,
                source_database text
            )""")

        con.execute("""\
            create table proteins (
                biogrid_id_interactor text,
                entrez_gene_interactor text,
                systematic_name_interactor text,
                official_symbol_interactor text,
                synonyms_interactor text,
                organism_interactor text
            )""")

        proteins = {}

        # Values that go in the links table
        link_indices = [0, 3, 4, 11, 12, 13, 14, 17, 18, 19, 20, 21, 22, 23]
        # Values that go in the proteins table
        interactor_a_indices = [3, 1, 5, 7, 9, 15]
        # Values that go in the proteins table
        interactor_b_indices = [4, 2, 6, 8, 10, 16]

        def to_none(val):
            return None if val == "-" else val

        def processlinks(rowiter):
            for row in rowiter:
                fields = map(to_none, row)
                yield [fields[i] for i in link_indices]

                interactor_a = [fields[i] for i in interactor_a_indices]
                interactor_b = [fields[i] for i in interactor_b_indices]
                proteins[interactor_a[0]] = interactor_a
                proteins[interactor_b[0]] = interactor_b

        con.executemany("""
            insert into links values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, processlinks(rows))

        con.executemany("""\
            insert into proteins values (?, ?, ?, ?, ?, ?)
            """, proteins.itervalues())
        con.commit()
        con.close()

    def init_db_index(self):
        """
        Will create an indexes (if not already present) in the database
        for faster searching by primary ids.

        """
        self.db.execute("""\
        create index if not exists index_on_biogrid_id_interactor_a
           on links (biogrid_id_interactor_a)
        """)
        self.db.execute("""\
        create index if not exists index_on_biogrid_id_interactor_b
           on links (biogrid_id_interactor_b)
        """)
        self.db.execute("""\
        create index if not exists index_on_biogrid_id_interactor
           on proteins (biogrid_id_interactor)
        """)


STRINGInteraction = namedtuple(
    "STRINGInteraciton",
    ["protein_id1",
     "protein_id2",
     "combined_score",
     "mode",
     "action",
     "score"]
)


class STRING(PPIDatabase):
    """
    Access `STRING <http://www.string-db.org/>`_ PPI database.
    """

    DATABASE_SCHEMA = """\
    Database schema
    ---------------
    table `links`:
        - `protein_id1`: id (text)
        - `protein_id2`: id (text)
        - `score`: combined score (int)

    table `actions`:
        - `protein_id1`: id (text)
        - `protein_id2`: id (text)
        - `mode`: mode (text)
        - `action`: action type (text)
        - `score`: action score (int)

    table `proteins`:
        - `protein_id`: protein id in STRING (text) (in the form of {taxid}.{name})
        - `taxid`: organism taxid (text)

    table `aliases`:
        - `protein_id: id (text)
        - `alias`: protein alias (text)
        - `source`: protein alias source (text)

    """
    DOMAIN = "PPI"
    FILENAME = "string.protein.{taxid}.sqlite"
    VERSION = "3.0"

    # Mapping from taxonomy.common_taxids() to taxids in STRING.

    TAXID_MAP = {"352472": "44689",  # Dictyostelium discoideum
                 "562": None,
                 "2104": "272634",   # Mycoplasma pneumoniae M129
                 "4530": "39947",    # Oryza sativa Japonica Group
                 "4754": None,
                 "8355": None,
                 "4577": None,
                 "11103": None}

    def __init__(self, taxid=None, database=None):
        if taxid is not None and database is not None:
            raise ValueError("taxid and database parameters are exclusive.")

        self.db = None

        if taxid is None and database is not None:
            if isinstance(database, sqlite3.Connection):
                self.db = database
                self.filename = None
            else:
                self.filename = database
                self.db = sqlite3.connect(database)
        elif taxid is not None and database is None:
            self.filename = serverfiles.localpath_download(
                self.DOMAIN, self.FILENAME.format(taxid=taxid)
            )
        elif taxid is None and database is None:
            # Back compatibility
            self.filename = serverfiles.localpath_download(
                "PPI", "string-protein.sqlite")
        else:
            assert False, "Not reachable"

        if self.db is None:
            self.db = sqlite3.connect(self.filename)

    @classmethod
    def default_db_filename(cls, taxid):
        return serverfiles.localpath(
            cls.DOMAIN, cls.FILENAME.format(taxid=taxid))

    @classmethod
    def common_taxids(cls):
        common = taxonomy.common_taxids()
        common = [cls.TAXID_MAP.get(id, id) for id in common]
        return filter(None, common)

    def organisms(self):
        """
        Return all organism taxids contained in this database.
        """
        cur = self.db.execute("select distinct taxid from proteins")
        return [r[0] for r in cur.fetchall()]

    def ids(self, taxid=None):
        """
        Return a list of all protein ids. If `taxid` is not None limit
        the results to ids from this organism only.

        """
        if taxid is not None:
            cur = self.db.execute("""\
                select protein_id
                from proteins
                where taxid=?
                """, (taxid,))
        else:
            cur = self.db.execute("""\
                select protein_id
                from proteins
                """)
        return [r[0] for r in cur.fetchall()]

    def synonyms(self, id):
        """
        Return a list of synonyms for primary `id` as reported by STRING
        (proteins.aliases.{version}.txt file)
        """
        cur = self.db.execute("""\
            select alias
            from aliases
            where protein_id=?
            """, (id,))
        res = cur.fetchall()
        return [r[0] for r in res]

    def synonyms_with_source(self, id):
        """
        Return a list of synonyms for primary `id` along with its
        source as reported by STRING (proteins.aliases.{version}.txt file)

        """
        cur = self.db.execute("""\
            select alias, source
            from aliases
            where protein_id=?
            """, (id,))
        res = cur.fetchall()
        return [(syn, set(source.split(" "))) \
                for syn, source in res]

    def all_edges(self, taxid=None):
        """
        Return a list of all edges. If taxid is not None return the edges
        for this organism only.

        .. note:: This may take some time (and memory).

        """
        if taxid is not None:
            cur = self.db.execute("""\
                select links.protein_id1, links.protein_id2, score
                from links join proteins on
                    links.protein_id1=proteins.protein_id
                where taxid=?
                """, (taxid,))
        else:
            cur = self.db.execute("""\
                select protein_id1, protein_id1, score
                from links
                """)
        return cur.fetchall()

    def edges(self, id):
        """
        Return a list of all edges (a list of 3-tuples (id1, id2, score)).
        """
        cur = self.db.execute("""\
            select protein_id1, protein_id2, score
            from links
            where protein_id1=?
            """, (id,))
        return cur.fetchall()

    def all_edges_annotated(self, taxid=None):
        res = []
        for id in self.ids(taxid):
            res.extend(self.edges_annotated(id))
        return res

    def edges_annotated(self, id):
        cur = self.db.execute("""\
            select links.protein_id1, links.protein_id2, links.score,
                   actions.action, actions.mode, actions.score
            from links left join actions on
                   links.protein_id1=actions.protein_id1 and
                   links.protein_id2=actions.protein_id2
            where links.protein_id1=?
        """, (id,))
        return map(STRINGInteraction._make, cur.fetchall())

    def search_id(self, name, taxid=None):
        if taxid is None:
            cur = self.db.execute("""\
                select proteins.protein_id
                from proteins natural join aliases
                where aliases.alias=?
            """, (name,))
        else:
            cur = self.db.execute("""\
                select proteins.protein_id
                from proteins natural join aliases
                where aliases.alias=? and proteins.taxid=?
            """, (name, taxid))
        return map(itemgetter(0), cur)

    @classmethod
    def download_data(cls, version, taxids=None):
        """
        Download the  PPI data for local work (this may take some time).
        Pass the version of the  STRING release e.g. v9.1.
        """
        if taxids is None:
            taxids = cls.common_taxids()

        for taxid in taxids:
            cls.init_db(version, taxid)

    @classmethod
    def init_db(cls, version, taxid, cache_dir=None, dbfilename=None):
        if cache_dir is None:
            cache_dir = serverfiles.localpath(cls.DOMAIN)

        if dbfilename is None:
            dbfilename = cls.default_db_filename(taxid)

        pjoin = os.path.join

        base_url = "http://string-db.org/newstring_download/"

        def paths(flatfile):
            url = "{flatfile}.{version}/{taxid}.{flatfile}.{version}.txt.gz"
            url = url.format(flatfile=flatfile, version=version, taxid=taxid)
            return posixpath.basename(url), base_url + url

        def ffname(pattern):
            return pattern.format(taxid=taxid, version=version)

        links_filename, links_url = paths("protein.links")

        actions_filename, actions_url = paths("protein.actions")

        aliases_filename, aliases_url = paths("protein.aliases")

        def download(filename, url):
            with open(pjoin(cache_dir, filename + ".tmp"), "wb") as dest:
                wget(url, dst_obj=dest, progress=True)

            shutil.move(pjoin(cache_dir, filename + ".tmp"),
                        pjoin(cache_dir, filename))

        for fname, url in [(links_filename, links_url),
                           (actions_filename, actions_url),
                           (aliases_filename, aliases_url)]:
            if not os.path.exists(pjoin(cache_dir, fname)):
                download(fname, url)

        links_fileobj = open(pjoin(cache_dir, links_filename), "rb")
        actions_fileobj = open(pjoin(cache_dir, actions_filename), "rb")
        aliases_fileobj = open(pjoin(cache_dir, aliases_filename), "rb")

        links_file = gzip.GzipFile(fileobj=links_fileobj)
        actions_file = gzip.GzipFile(fileobj=actions_fileobj)
        aliases_file = gzip.GzipFile(fileobj=aliases_fileobj)

        progress = ConsoleProgressBar("Processing {}:".format(links_filename))
        progress(0.0)

        def st_size(filename):
            return os.stat(pjoin(cache_dir, filename)).st_size

        filesize = st_size(links_filename)

        con = sqlite3.connect(dbfilename)

        with con:
            cls.clear_db(con)

            links_file.readline()  # read the header line

            reader = csv.reader(links_file, delimiter=" ")

            def read_links(reader, progress):
                for i, (p1, p2, score) in enumerate(reader):
                    yield p1, p2, int(score)

                    if i % 100000 == 0:
                        # Update the progress every 100000 lines
                        progress(100.0 * links_fileobj.tell() / filesize)

            con.executemany("INSERT INTO links VALUES (?, ?, ?)",
                            read_links(reader, progress))

            progress.finish()

            def part(string, sep, part):
                return string.split(sep)[part]

            con.create_function("part", 3, part)
            con.execute("""
                INSERT INTO proteins
                SELECT protein_id1, part(protein_id1, '.', 0)
                FROM (SELECT DISTINCT(protein_id1)
                     FROM links
                     ORDER BY protein_id1)
            """)

            filesize = st_size(actions_filename)

            actions_file.readline()  # read header line

            progress = ConsoleProgressBar("Processing actions:")
            reader = csv.reader(actions_file, delimiter="\t")

            def read_actions(reader):
                for i, (p1, p2, mode, action, a_is_acting, score) in \
                        enumerate(reader):
                    yield p1, p2, mode, action, int(score)

                    if i % 10000 == 0:
                        progress(100.0 * actions_fileobj.tell() / filesize)

            con.executemany("INSERT INTO actions VALUES (?, ?, ?, ?, ?)",
                            read_actions(reader))

            progress.finish()

            filesize = st_size(aliases_filename)
            aliases_file.readline()  # read header line

            progress = ConsoleProgressBar("Processing aliases:")

            reader = csv.reader(aliases_file, delimiter="\t")

            def read_aliases(reader, progress):
                for i, (taxid, name, alias, source) in enumerate(reader):
                    yield (".".join([taxid, name]),
                           alias.decode("utf-8", errors="ignore"),
                           source.decode("utf-8", errors="ignore"))
                    if i % 10000 == 0:
                        progress(100.0 * aliases_fileobj.tell() / filesize)

            con.executemany("INSERT INTO aliases VALUES (?, ?, ?)",
                            read_aliases(reader, progress))

            progress.finish()

            print "Indexing the database"
            cls.create_db_index(con)

            con.executescript("""
                DROP TABLE IF EXISTS version;
                CREATE TABLE version (
                     string_version text,
                     api_version text
                );""")

            con.execute("""
                INSERT INTO version
                VALUES (?, ?)""", (version, cls.VERSION))

    @classmethod
    def clear_db(cls, dbcon):
        dbcon.executescript(textwrap.dedent("""
            DROP TABLE IF EXISTS links;
            DROP TABLE IF EXISTS proteins;
            DROP TABLE IF EXISTS actions;
            DROP TABLE IF EXISTS aliases;
        """))

        dbcon.executescript(textwrap.dedent("""
            CREATE TABLE links
                (protein_id1 TEXT, protein_id2 TEXT, score INT);

            CREATE TABLE proteins
                (protein_id TEXT, taxid TEXT);

            CREATE TABLE actions
                (protein_id1 TEXT, protein_id2 TEXT, mode TEXT,
                 action TEXT, score INT);

            CREATE TABLE aliases
                (protein_id TEXT, alias TEXT, source TEXT);
        """))

    @classmethod
    def create_db_index(cls, dbcon):
        dbcon.executescript(textwrap.dedent("""
            CREATE INDEX IF NOT EXISTS index_link_protein_id1
                ON links (protein_id1);

            CREATE INDEX IF NOT EXISTS index_action_protein_id1
                ON actions (protein_id1);

            CREATE INDEX IF NOT EXISTS index_proteins_id
                ON proteins (protein_id);

            CREATE INDEX IF NOT EXISTS index_taxids
                ON proteins (taxid);

            CREATE INDEX IF NOT EXISTS index_aliases_id
                ON aliases (protein_id);

            CREATE INDEX IF NOT EXISTS index_aliases_alias
                ON aliases (alias);
        """))


STRINGDetailedInteraction = namedtuple(
    "STRINGDetailedInteraction",
    ["protein_id1",
     "protein_id2",
     "combined_score",
     "mode",
     "action",
     "score",
     "neighborhood",
     "fusion",
     "cooccurence",
     "coexpression",
     "experimental",
     "database",
     "textmining"]
)


class STRINGDetailed(STRING):
    """
    Access `STRING <http://www.string-db.org/>`_ PPI database.
    This class also allows access to subscores per channel.

    .. note::
        This data is released under a `Creative Commons
        Attribution-Noncommercial-Share Alike 3.0 License
        <http://creativecommons.org/licenses/by-nc-sa/3.0/>`_.

        If you want to use this data for commercial purposes you must
        get a license from STRING.

    """

    DATABASE_SCHEMA = """\
    DATABASE SCHEMA
    ===============

    table `evidence`:
        - `protein_id1`: protein id (text)
        - `protein_id2`: protein id (text)
        - `neighborhood`: score (int)
        - `fusion`: score (int)
        - `cooccurence`: score (int)
        - `coexpression`: score (int)
        - `experimental`: score (int)
        - `database`: score (int)
        - `textmining`: score (int)

    """
    FILENAME_DETAILED = "string.protein.detailed.{taxid}.sqlite"

    def __init__(self, taxid=None, database=None, detailed_database=None):
        STRING.__init__(self, taxid, database)
        if taxid is not None and detailed_database is not None:
            raise ValueError("taxid and detailed_database are exclusive")

        db_file = serverfiles.localpath(self.DOMAIN, self.FILENAME)
        if taxid is not None and detailed_database is None:
            detailed_database = serverfiles.localpath_download(
                self.DOMAIN,
                self.FILENAME_DETAILED.format(taxid=taxid)
            )
        elif taxid is None and detailed_database is not None:
            detailed_database = detailed_database
        elif taxid is None and detailed_database is None:
            # Back compatibility
            detailed_database = serverfiles.localpath_download(
                "PPI", "string-protein-detailed.sqlite")

        self.db_detailed = sqlite3.connect(detailed_database)
        self.db_detailed.execute("ATTACH DATABASE ? as string", (db_file,))

    def edges_annotated(self, id):
        edges = STRING.edges_annotated(self, id)
        edges_nc = []
        for edge in edges:
            id1, id2 = edge.protein_id1, edge.protein_id2
            cur = self.db_detailed.execute("""
                SELECT neighborhood, fusion, cooccurence, coexpression,
                       experimental, database, textmining
                FROM evidence
                WHERE protein_id1=? AND protein_id2=?
                """, (id1, id2))
            res = cur.fetchone()
            if res:
                evidence = res
            else:
                evidence = [0] * 7
            edges_nc.append(
                STRINGDetailedInteraction(*(tuple(edge) + tuple(evidence)))
            )
        return edges_nc

    @classmethod
    def init_db(cls, version, taxid, cache_dir=None, dbfilename=None):
        if cache_dir is None:
            cache_dir = serverfiles.localpath(cls.DOMAIN)
        if dbfilename is None:
            dbfilename = serverfiles.localpath(
                cls.DOMAIN,
                "string-protein-detailed.{taxid}.sqlite".format(taxid=taxid)
            )

        pjoin = os.path.join

        base_url = "http://string-db.org/newstring_download/"
        filename = "{taxid}.protein.links.detailed.{version}.txt.gz"
        filename = filename.format(version=version, taxid=taxid)
        url = base_url + "protein.links.detailed.{version}/" + filename
        url = url.format(version=version)

        if not os.path.exists(pjoin(cache_dir, filename)):
            wget(url, cache_dir, progress=True)

        links_fileobj = open(pjoin(cache_dir, filename), "rb")
        links_file = gzip.GzipFile(fileobj=links_fileobj)

        con = sqlite3.connect(dbfilename)
        with con:
            con.execute("""
                DROP TABLE IF EXISTS evidence
            """)

            con.execute("""
                CREATE TABLE evidence(
                     protein_id1 TEXT,
                     protein_id2 TEXT,
                     neighborhood INTEGER,
                     fusion INTEGER,
                     cooccurence INTEGER,
                     coexpression INTEGER,
                     experimental INTEGER,
                     database INTEGER,
                     textmining INTEGER
                    )
                """)

            links = csv.reader(links_file, delimiter=" ")
            links.next()  # Read header
            filesize = os.stat(pjoin(cache_dir, filename)).st_size

            progress = ConsoleProgressBar("Processing links file:")
            progress(1.0)

            def read_links(reader):
                for i, (p1, p2, n, f, c, cx, ex, db, t, _) in \
                        enumerate(reader):
                    yield p1, p2, n, f, c, cx, ex, db, t

                    if i % 10000 == 0:
                        progress(100.0 * links_fileobj.tell() / filesize)

            con.executemany("""
                INSERT INTO evidence
                VALUES  (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, read_links(links))

            progress.finish()

            print "Indexing"
            con.execute("""\
                CREATE INDEX IF NOT EXISTS index_evidence
                    ON evidence (protein_id1, protein_id2)
            """)

            con.executescript("""
                DROP TABLE IF EXISTS version;

                CREATE TABLE version (
                     string_version text,
                     api_version text
                );
                """)

            con.execute("""
                INSERT INTO version
                VALUES (?, ?)""", (version, cls.VERSION))

    @classmethod
    def download_data(cls, version, taxids=None):
        if taxids is None:
            taxids = cls.common_taxids()

        for taxid in taxids:
            cls.init_db(version, taxid)


##########
# Obsolete
##########

class Interaction(object):
    def __init__(self, protein1, protein2, ref1=None, ref2=None, conf1=None, conf2=None):
        self.protein1, self.protein2 = protein1, protein2
        self.ref1, self.ref2 = ref1, ref2
        self.conf1, self.conf2 = conf1, conf2
        self.org1, self.org2 = None, None


class MIPS(object):
    VERSION = 1

    def __init__(self):
        self.load()

    def load(self):
        self.protein_names = defaultdict(set)
        self.refs = {}
        self.confidance = {}

        def process(element):
            d = {}
            participants = element.getElementsByTagName("proteinParticipant")
            proteins = []
            for protein in participants:
                interactor = protein.getElementsByTagName("proteinInteractor")[0]
                names = []
                for name in interactor.getElementsByTagName("shortLabel") + \
                            interactor.getElementsByTagName("fullName"):
                    names.append((name.tagName, name.childNodes[0].data))

                refs = []
                for ref in interactor.getElementsByTagName("primaryRef"):
                    refs += [(ref.tagName, ref.attributes.items())]
                org = dict(interactor.getElementsByTagName("organism")[0].attributes.items()).get("ncbiTaxId")
                conf = protein.getElementsByTagName("confidence")[0].attributes.items()
                proteins.append((names, refs, conf, org))
            interaction = Interaction(proteins[0][0][1][1], proteins[1][0][1][1])
            interaction.ref1, interaction.ref2 = proteins[0][1], proteins[1][1]
            interaction.conf1, interaction.conf2 = proteins[0][2], proteins[1][2]
            interaction.org1, interaction.org2 = proteins[0][3], proteins[1][3]

            self.protein_names[interaction.protein1].add(proteins[0][0][0][1])
            self.protein_names[interaction.protein2].add(proteins[1][0][0][1])

            return interaction

        document = minidom.parse(serverfiles.localpath_download("PPI", "allppis.xml"))
        interactions = document.getElementsByTagName("interaction")
        self.interactions = [process(interaction) for interaction in interactions]

        self.protein_interactions = defaultdict(set)

        for inter in self.interactions:
            self.protein_names[inter.protein1] = dict(inter.ref1[0][1]).get("id")
            self.protein_names[inter.protein2] = dict(inter.ref2[0][1]).get("id")
            self.protein_interactions[inter.protein1].add(inter)
            self.protein_interactions[inter.protein2].add(inter)

    def __iter__(self):
        return iter(self.interactions)

    @classmethod
    def download(cls):
        src = urllib2.urlopen("http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz")
        dest = serverfiles.localpath("PPI", "mppi.gz")
        shutil.copyfileobj(src, open(dest, "wb"))

    @classmethod
    @pickled_cache(None, [("PPI", "allppis.xml")], version=1)
    def _get_instance(cls):
        return MIPS()

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_instance"):
            cls._instance = cls._get_instance()
        return cls._instance


def mips_interactions(protein=None):
    mips = MIPS.get_instance()
    if protein is None:
        return list(mips)
    else:
        return mips.protein_interactions.get(protein)


def mips_proteins():
    return set(MIPS.get_instance().protein_names.keys())


class BioGRIDInteraction(object):
    """ An object representing a BioGRID interaction. Each member of this object
    represents a data from a single column of BIOGRID-ALL.tab file.

    Attributes:
        - *interactor_a*    - BioGRID identifier
        - *interactor_b*    - BioGRID identifier
        - *official_symbol_a*    - An official symbol for *interactor_a*
        - *official_symbol_b*    - An official symbol for *interactor_b*
        - *aliases_for_a*    - Aliases separated by '|'
        - *aliases_for_b*    - Aliases separated by '|'
        - *experimental_system*     - Experimental system (see BioGRID documentation on www.thebiogrid.org for a list of valid entrys)
        - *source*    -
        - *organism_a_id*    - NCBI Taxonomy identifier for *interactor_a*'s organism
        - *organism_b_id*    - NCBI Taxonomy identifier for *interactor_b*'s organism
    """
    __slots__ = ["interactor_a", "interactor_b", "official_symbol_a", "official_symbol_b", "aliases_for_a", "aliases_for_b", "experimental_system", "source", "pubmed_id", "organism_a_id", "organism_b_id"]

    def __init__(self, line):
        for attr, val in zip(self.__slots__, line.split("\t")):
            setattr(self, attr, val)


class _BioGRID_Old(object):
    """ A BioGRID database interface
    Example::
        >>> ## finding all interactions for Homo sapiens sapiens
        >>> grid = BioGRID(case_insensitive=True)
        >>> proteins = proteins = biogrid.proteins() ## All proteins
        >>> proteins = [p for p in proteins if any(["9606" in [int.organism_a_id, int.organism_b_id] for int in grid.get(p)])]
    """
    VERSION = 1

    def __init__(self, case_insensitive=True):
#        warnings.warn("obiPPi._BioGRID_Old class is deprecated. Use obiPPI.BioGRID")
        self.case_insensitive = case_insensitive
        self._case = (lambda name: name.lower()) if self.case_insensitive else (lambda name: name)
        self.load()

    def load(self):
        text = open(serverfiles.localpath_download("PPI", "BIOGRID-ALL.tab"), "rb").read()
        text = text.split("SOURCE\tPUBMED_ID\tORGANISM_A_ID\tORGANISM_B_ID\n", 1)[-1]
        self.interactions = [BioGRIDInteraction(line) for line in text.splitlines() if line.strip()]

        self.protein_interactions = defaultdict(set)
        self.protein_names = {}

        case = self._case

        def update(keys, value, collection):
            for k in keys:
                collection.setdefault(k, set()).add(value)

        for inter in self.interactions:
            update(map(case, [inter.official_symbol_a] + inter.aliases_for_a.split("|")), case(inter.interactor_a), self.protein_names)
            update(map(case, [inter.official_symbol_b] + inter.aliases_for_b.split("|")), case(inter.interactor_b), self.protein_names)

            self.protein_interactions[case(inter.interactor_a)].add(inter)
            self.protein_interactions[case(inter.interactor_b)].add(inter)

        self.protein_interactions = dict(self.protein_interactions)

        if case("N/A") in self.protein_names:
            del self.protein_names[case("N/A")]

    def proteins(self):
        """ Return all protein names in BioGRID (from INTERACTOR_A, and INTERACTOR_B columns).
        """
        return self.protein_interactions.keys()

    def __iter__(self):
        """ Iterate over all BioGRIDInteraction objects.
        """
        return iter(self.interactions)

    def __getitem__(self, key):
        """ Return a list of protein interactions that a protein is a part of.
        """
        key = self._case(key)
#        keys = self.protein_alias_matcher.match(key)
        if key not in self.protein_interactions:
            keys = self.protein_names.get(key, [])
        else:
            keys = [key]
        if keys:
            return list(reduce(set.union, [self.protein_interactions.get(k, []) for k in keys], set()))
        else:
            raise KeyError(key)

    def get(self, key, default=None):
        """ Return a list of protein interactions that a protein is a part of.
        """
        key = self._case(key)
#        keys = self.protein_alias_matcher.match(key)
        if key not in self.protein_interactions:
            keys = self.protein_names.get(keys, [])
        else:
            keys = [key]
        if keys:
            return list(reduce(set.union, [self.protein_interactions.get(k, []) for k in keys], set()))
        else:
            return default

    @classmethod
    def get_instance(cls):
        if getattr(cls, "_instance", None) is None:
            cls._instance = _BioGRID_Old()
        return cls._instance


def biogrid_interactions(name=None):
    """Return a list of protein interactions (BioGRIDInteraction objects) that a protein is a part of.
    """
    if name:
        return list(_BioGRID_Old.get_instance().get(name, set()))
    else:
        return _BioGRID_Old.get_instance().interactions


def biogrid_proteins():
    """ Return all protein names in BioGRID (from INTERACTOR_A, and INTERACTOR_B columns).
    """
    return _BioGRID_Old.get_instance().proteins()


if __name__ == "__main__":
    for protein in mips_proteins():
        print "Protein", protein, "interacts with",
        print ",".join(set(reduce(list.__add__, [[inter.protein1, inter.protein2] for inter in mips_interactions(protein)], [])) - set([protein]))

