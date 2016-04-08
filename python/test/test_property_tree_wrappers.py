# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from pycap import PropertyTree
from os import remove
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

    def test_get_children(self):
        ptree = PropertyTree()
        # put child
        child = PropertyTree()
        child.put_double('prune', 6.10)
        ptree.put_child('a.g', child)
        self.assertEqual(ptree.get_double('a.g.prune'), 6.10)
        # get child
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

    def test_parsers(self):
        # populate a populate property tree from an input file
        ptree = PropertyTree()
        # INFO parser
        info_file = 'input.info'
        today = 'april 8th, 2016'
        with open(info_file, 'w') as fout:
            fout.write('date "'+today+'"')
        ptree.parse_info(info_file)
        self.assertEqual(ptree.get_string('date'), today)
        remove(info_file)
        # XML parser
        xml_file = 'input.xml'
        pi = 3.14
        with open(xml_file, 'w') as fout:
            fout.write('<pi>'+str(pi)+'</pi>')
        ptree.parse_xml(xml_file)
        self.assertEqual(ptree.get_double('pi'), pi)
        remove(xml_file)
        # JSON parser
        json_file = 'input.json'
        with open(json_file, 'w') as fout:
            fout.write('{')
            fout.write('  "foo":')
            fout.write('  {')
            fout.write('    "bar": false')
            fout.write('  }')
            fout.write('}')
        ptree.parse_json(json_file)
        self.assertFalse(ptree.get_bool('foo.bar'))
        remove(json_file)

if __name__ == '__main__':
    unittest.main()
