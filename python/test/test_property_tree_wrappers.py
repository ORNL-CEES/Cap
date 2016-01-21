from pycap import PropertyTree
import unittest


class boostPropertyTreePythonWrappersTestCase(unittest.TestCase):
    def test_pickle_support(self):
        import pickle
        src = PropertyTree()
        src.put_double('pi', 3.14)
        src.put_string('greet', 'bonjour')
        p = pickle.dumps(src)
        dst = pickle.loads(p)
        self.assertEqual(dst.get_double('pi'), 3.14)
        self.assertEqual(dst.get_string('greet'), 'bonjour')

    def test_get_with_default_value(self):
        ptree = PropertyTree()
        # double
        self.assertEqual(ptree.get_double_with_default_value('missing_double',
                                                             3.14), 3.14)
        ptree.put_double('present_double', 1.41)
        self.assertEqual(ptree.get_double_with_default_value('present_double',
                                                             3.14), 1.41)
        # string
        self.assertEqual(ptree.get_string_with_default_value('missing_string',
                                                             'missing'),
                         'missing')
        ptree.put_string('present_string', 'present')
        self.assertEqual(ptree.get_string_with_default_value('present_string',
                                                             'missing'),
                         'present')
        # int
        self.assertEqual(ptree.get_int_with_default_value('missing_int', 255),
                         255)
        ptree.put_int('present_int', 0)
        self.assertEqual(ptree.get_int_with_default_value('present_int', 255),
                         0)
        # bool
        self.assertEqual(ptree.get_bool_with_default_value('missing_bool',
                                                           True), True)
        ptree.put_bool('present_bool', False)
        self.assertEqual(ptree.get_bool_with_default_value('present_bool',
                                                           True), False)

    def test_get_array(self):
        ptree = PropertyTree()
        # array of double
        ptree.put_string('array_double', '3.14,1.41')
        array_double = ptree.get_array_double('array_double')
        self.assertEqual(array_double, [3.14, 1.41])
        # ... string
        ptree.put_string('array_int', '1,2,3')
        array_int = ptree.get_array_int('array_int')
        self.assertEqual(array_int, [1, 2, 3])
        # ... int
        ptree.put_string('array_string', 'uno,dos,tres,cuatro')
        array_string = ptree.get_array_string('array_string')
        self.assertEqual(array_string, ['uno', 'dos', 'tres', 'cuatro'])
        # ... bool
        ptree.put_string('array_bool', 'true,FALSE,False')
        array_bool = ptree.get_array_bool('array_bool')
        self.assertEqual(array_bool, [True, False, False])

    def test_property_tree(self):
        # ptree as container to store int, double, string, and bool
        ptree = PropertyTree()
        ptree.put_int('dim', 3)
        self.assertEqual(ptree.get_int('dim'), 3)
        ptree.put_double('path.to.pi', 3.14)
        self.assertEqual(ptree.get_double('path.to.pi'), 3.14)
        ptree.put_string('good.news', 'it works')
        self.assertEqual(ptree.get_string('good.news'), 'it works')
        ptree.put_bool('is.that.a.good.idea', False)
        self.assertEqual(ptree.get_bool('is.that.a.good.idea'), False)

    def test_get_child(self):
        ptree = PropertyTree()
        ptree.put_string('child.name', 'clement')
        ptree.put_int('child.age', -2)
        child = ptree.get_child('child')
        self.assertEqual(child.get_string('name'), 'clement')
        self.assertEqual(child.get_int('age'), -2)

    def test_raise_exceptions(self):
        ptree = PropertyTree()
        # property tree will throw if the specified path does not exist
        self.assertRaises(RuntimeError, ptree.get_int, 'path.does.not.exist')
        # or if the translation fails
        ptree.put_string('some.path.to.a.string', 'not a double')
        self.assertRaises(RuntimeError, ptree.get_double,
                          'some.path.to.a.string')

    # TODO
    def test_parse(self):
        ptree = PropertyTree()
        ptree.parse_xml('device.xml')

if __name__ == '__main__':
    unittest.main()
