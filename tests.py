"""nosetests for project"""
import nose.tools
import utils, random

def test_w_stat():
     """nosetest for utils.w_stat"""

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

def test_my_wil_rank_sum_gtr():
     """nosetest for utils.my_wil_rank_sum_gtr"""

     # test gtr dist
     vals = []
     for x in xrange(1000):
          vals.append(random.randint(0, 10000))
     vals.sort()
     low_vals = vals[0:300]
     high_vals = vals[800:]

     pval = utils.my_wil_rank_sum_gtr(high_vals, low_vals,
                                      vals)
     nose.tools.assert_true(pval < float(.01))

     # test eq dist
     low_vals = []
     for x in xrange(300):
          low_vals.append( vals[random.randint(0,999)] )
     high_vals = []
     for x in xrange(300):
          high_vals.append( vals[random.randint(0,999)] )

     pval = utils.my_wil_rank_sum_gtr(high_vals, low_vals,
                                      vals)
     nose.tools.assert_true(pval > float(.05),
                            'pval was significant, but should not be '
                            + str(pval))
