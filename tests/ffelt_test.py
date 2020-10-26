import unittest
from binext import ffelt

class testAdd (unittest.TestCase):

    def test_sum(self):
        sums16={
            (None, None): None, (1, 1): None, (13, 13): None, (15, 15): None,
            (1, None): 1, (None, 1): 1,
            (13, None): 13, (None, 13): 13,
            (15, None): 0, (None, 15): 0
        }

        for a,b in sums16:
            gfa=ffelt(a,16)
            gfb=ffelt(b,16)
            gfc=ffelt(sums16[(a,b)],16)
            self.assertEqual(gfa + gfb, gfc,
                             str(gfa)+"+"+str(gfb)+"should be"+str(gfc)
                             )

if __name__ == '__main__':
    unittest.main()
