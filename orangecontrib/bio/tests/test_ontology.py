import doctest
import unittest

from six import StringIO

from orangecontrib.bio import ontology

class TestOntology(unittest.TestCase):
    def test_oboobject(self):
        stanza = '''[Term]
id: FOO:001
name: bar
'''
        term = ontology.OBOObject.parse_stanza(stanza)

#        self.assertEqual(term.get_value("id"), "FOO:001")
#        self.assertEqual(term.get_value("name"), "bar")

        term.add_tag("tt", "3")

    def test_ontology(self):
        seinfeld = StringIO("""
[Typedef]
id: parent

[Typedef]
id: child
inverse_of: parent ! not actually used yet

[Term]
id: 001
name: George

[Term]
id: 002
name: Estelle
relationship: parent 001 ! George

[Term]
id: 003
name: Frank
relationship: parent 001 ! George

""")  # TODO: fill the ontology with all characters
        stanza = '''[Term]
id: FOO:001
name: bar
'''
#        term = ontology.OBOObject.parse_stanza(stanza)

        seinfeld = ontology.OBOOntology(seinfeld)
#        print(seinfeld.child_edges("001"))

def load_tests(loader, tests, ignore):
    stanza = '''[Term]
id: FOO:001
name: bar
'''
    term = ontology.OBOObject.parse_stanza(stanza)
    tests.addTests(doctest.DocTestSuite(ontology,
                                        extraglobs={"term": term},
                                        optionflags=doctest.ELLIPSIS))
    return tests
