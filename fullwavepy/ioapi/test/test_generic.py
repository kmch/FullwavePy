import numpy as np
from unittest import TestCase, skip
from fullwavepy.ioapi.generic import read_dict

class Test(TestCase):
  def test_read_dict(self):
    fname = 'test-Dict.key'
    with open(fname, 'w') as f:
      f.write('a : 1\n')
      f.write('\n') # empty line
      f.write('b b : 2\n')
      f.write('akfjkajf\n') # bad line
      f.write('!c : 3\n')
      f.write('MAX TIME : 1000')
    d = read_dict(fname)
    # print(d)
    assert d == {'a': '1', 'b b': '2', 'ttime': '1.0'}