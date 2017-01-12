import sys
import unittest

import Orange


@unittest.skipIf(sys.version_info < (3,), "Python 3 only")
class TestGenapi(unittest.TestCase):
    def test_login(self):
        from orangecontrib.bio import resolwe
        email = 'anonymous@genialis.com'
        password = 'anonymous'
        url = 'https://dictyexpress.research.bcm.edu'

        self.assertTrue(resolwe.connect(email, password, url, 'genesis'))
        self.assertRaises(Exception, resolwe.connect, email, "123", url, 'genesis')
        self.assertRaises(Exception, resolwe.connect, email, password, 'testUrl', 'genesis')

    def test_objects(self):
        from orangecontrib.bio import resolwe
        gen = resolwe.connect('anonymous@genialis.com', 'anonymous',
                              'https://dictyexpress.research.bcm.edu', 'genesis')
        etc_objects = gen.fetch_etc_objects()

        self.assertEqual(type(etc_objects), list)  # test if return type is list
        self.assertTrue(etc_objects)  # test if list is not empty
        self.assertEqual(etc_objects[0].type, 'data:etc:')  # test if it contains correct objects

        #  Test experiment D. purpureu

        for obj in etc_objects:
            if obj.name == 'D. purpureum':
                test_experiment = obj

        self.assertTrue(test_experiment)
        self.assertEqual(test_experiment.id, '564a509e6b13390ffb40d4c8')

        json = gen.download_etc_data(test_experiment.id)

        self.assertEqual(type(json), dict)
        self.assertEqual(len(json["etc"].keys()), 2)
        self.assertEqual(len(json["etc"]["genes"].keys()), 12410)
        self.assertEqual(json["etc"]["genes"]["DPU_G0071544"], [0.0, 0.1787337345055, 4.20485453935, 20.002575156149998,
                          19.52080354305, 18.7919080288, 12.38709403699])
        self.assertEqual(json["etc"]["timePoints"], [0, 4, 8, 12, 16, 20, 24])

        self.assertEqual(type(gen.etc_to_table(json)), Orange.data.table.Table)


if __name__ == '__main__':
    unittest.main()
