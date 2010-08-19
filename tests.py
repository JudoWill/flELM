"""nosetests for project"""
import nose.tools
import utils, random

def test_w_stat():
     """nosetest for w_stat"""

     # sorted vals
     ls1 = (1,2,3)
     ls2 = (4,5,6)
     w1, w2 = utils.w_stat(ls1, ls2)
     
     nose.tools.assert_equal(w1, 6)
     nose.tools.assert_equal(w2, 15)
     nose.tools.assert_equal(w1+w2, 21)

     # test sorting
     ls1 = (1,2,4)
     ls2 = (3,5,6)
     w1, w2 = utils.w_stat(ls1, ls2)
     
     nose.tools.assert_equal(w1, 7)
     nose.tools.assert_equal(w2, 14)
     nose.tools.assert_equal(w1+w2, 21)

     # test duplicates
     ls1 = (1,2,4)
     ls2 = (3,5,6,4)
     w1, w2 = utils.w_stat(ls1, ls2)

     nose.tools.assert_equal(w1, 7.5)
     nose.tools.assert_equal(w2, 20.5)
     nose.tools.assert_equal(w1+w2, 28)
          
     # test duplicates in same list
     ls1 = (1,2,4,4)
     ls2 = (3,5,6,4)
     w1, w2 = utils.w_stat(ls1, ls2)
     
     nose.tools.assert_equal(w1, 13)
     nose.tools.assert_equal(w2, 23)
     nose.tools.assert_equal(w1+w2, 36)

     
     
     
