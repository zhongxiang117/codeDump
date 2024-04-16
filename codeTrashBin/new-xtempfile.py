#!/usr/bin/env python3

import os
import tempfile


# this is a continuous testing for file xtempfile.py
# its link: https://bugs.python.org/file50004/xtempfile.py
# issue link: https://bugs.python.org/issue44013


def test_7():
    tmp1 = tempfile.TemporaryFile()
    os.close(tmp1.fileno())
    tmp2 = tempfile.TemporaryFile()
    del tmp1

    tmp2.seek(0)    # OSError: Bad file descriptor
                    # how can del tmp1 influence different file descriptor tmp2?

    # by checking source codes:
    #   https://github.com/python/cpython/blob/3.9/Lib/tempfile.py#L440
    #
    # del method is actually implemented by self.close() method
    #
    # Does the problem mean tempfile.TemporaryFile is a singleton class??
    # In other words, del method only closes the newest file descriptor.
    #
    # no, apparently it is not.



# please compare this with test_7
# only the order of "del" is different
def test_8():
    tmp1 = tempfile.TemporaryFile()
    os.close(tmp1.fileno())
    del tmp1        # technically, tmp1 is closed twice

    tmp2 = tempfile.TemporaryFile()
    tmp2.seek(0)    # works well



# please compare this with test_8
# del method is put to slightly below
def test_9():
    tmp1 = tempfile.TemporaryFile()
    os.close(tmp1.fileno())
    
    tmp2 = tempfile.TemporaryFile()
    tmp2.seek(0)    # works well

    del tmp1

    tmp2.write(b'good')     # works well

    # however, file descriptor tmp2 cannot do seek or flush or many operations
    # this is exactly the similar problem that happens in test_7

    #tmp2.flush
    #tmp2.seek(0)
    #tmp2.readline()
    tmp2.readlines()



# please carefully compare this with test_9
# os.close is commented out
# everything again works fine now
def test_10():
    tmp1 = tempfile.TemporaryFile()
    #os.close(tmp1.fileno())
    
    tmp2 = tempfile.TemporaryFile()
    tmp2.seek(0)    # works well

    del tmp1

    tmp2.write(b'good')     # works well

    # works well

    #tmp2.flush
    tmp2.seek(0)
    #tmp2.readline()
    print(tmp2.readlines())



# conclusion
# it should be a bug

# the problem seems like,
#
# for del method, if its original file descriptor is closed, e.g., by os.close()
# ideally, it should do nothing; however, from the testing result, it will
# search the just next file descriptor, and then close it if find.
#
# to test this assumption, three file descriptors are created.


# this is just like the assumed!
def test_11():
    tmp1 = tempfile.TemporaryFile()
    os.close(tmp1.fileno())     # tmp1 is closed
    
    tmp2 = tempfile.TemporaryFile()
    tmp2.seek(0)    # works well

    tmp3 = tempfile.TemporaryFile()
    tmp3.seek(0)    # works well

    del tmp1        # by assumption, this will close tmp2, but tmp3 should be OK


    tmp2.write(b'good')     # works well

    # file descriptor tmp2 cannot do any operations
    #tmp2.flush
    #tmp2.seek(0)
    #tmp2.readline()
    #print(tmp2.readlines())


    # file descriptor tmp3 is not influenced at all
    tmp3.write(b'good')     # works well
    
    #tmp3.flush
    tmp3.seek(0)
    #tmp3.readline()
    print(tmp3.readlines())



# possible fixes
#
# make __del__ check whether input file descriptor is closed or not
#
# if it is true: do nothing
# otherwise: close it
#
# However,
# it still cannot explain the problem the I had in test_4 (in file: xtempfile)



