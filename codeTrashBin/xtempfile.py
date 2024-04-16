import tempfile


# works well
def test_1():
    print('Note: inside test_1')
    fd = tempfile.TemporaryFile()
    print('test_1 =>',id(fd),id(fd.name))
    fd.write(b'test_1: my strings')
    fd.seek(0)
    with open(fd.name) as f:
        print('test_1: fd is opened and closed in with')
        print(f.readlines())
    
    # fd is already closed inside open
    # it is expected
    #
    #fd.close()


# same contents like test_1
# works well
def test_2():
    print('Note: inside test_2')
    fd = tempfile.TemporaryFile()
    print('test_2 =>',id(fd),id(fd.name))
    fd.write(b'test_2: my strings')
    fd.seek(0)
    with open(fd.name) as f:
        print('test_2: fd is opened and closed in with')
        print(f.readlines())



# functions combined test_1 & test_2
# works well
def test_3():
    print('Note: inside test_3')
    fd = tempfile.TemporaryFile()
    print('test_3 =>',id(fd),id(fd.name))
    fd.close()

    test_1()
    test_2()



# codes combined test_1 & test_2
# has problem
def test_4():
    print('Note: inside test_4')
    fd = tempfile.TemporaryFile()
    print('test_1 =>',id(fd),id(fd.name))
    fd.write(b'test_1: my strings')
    fd.seek(0)
    fd.flush()
    with open(fd.name) as f:
        print('test_1: fd is opened and closed in with')
        print(f.readlines())


    fd = tempfile.TemporaryFile()
    print('test_2 =>',id(fd),id(fd.name))
    fd.write(b'test_2: my strings')

    # following will have:
    #
    #   OSError: [Errno 9] Bad file descriptor
    #
    # for any types of operations, e.gs;
    # fd.close(),  fd.flush()
    #
    # the name fd cannot be reused even in new instance
    fd.seek(0)
    with open(fd.name) as f:
        print('test_2: fd is opened and closed in with')
        print(f.readlines())



# works well by reset file descriptor
def test_5():
    print('Note: inside test_5')
    fd = tempfile.TemporaryFile()
    print('test_1 =>',id(fd),id(fd.name))
    fd.write(b'test_1: my strings')
    fd.seek(0)
    fd.flush()
    with open(fd.name) as f:
        print('test_1: fd is opened and closed in with')
        print(f.readlines())
    

    # reset to fix
    fd = 'some thing else'

    fd = tempfile.TemporaryFile()
    print('test_2 =>',id(fd),id(fd.name))
    fd.write(b'test_2: my strings')
    fd.seek(0)
    with open(fd.name) as f:
        print('test_2: fd is opened and closed in with')
        print(f.readlines())



# works well by the help of with
def test_6():
    print('Note: inside test_6')
    fd = tempfile.TemporaryFile()
    print('test_1 =>',id(fd),id(fd.name))
    fd.write(b'test_1: my strings')
    fd.seek(0)
    fd.flush()
    with open(fd.name) as f:
        print('test_1: fd is opened and closed in with')
        print(f.readlines())


    # second with
    with tempfile.TemporaryFile() as fd:
        fd = tempfile.TemporaryFile()
        print('test_2 =>',id(fd),id(fd.name))
        fd.write(b'test_2: my strings')
        fd.seek(0)
        fd.read()

    # fd is closed in second with
    # it is expected
    #
    #fd.close()



# codes combined test_1 & test_2
# has problem
def test_4():
    print('Note: inside test_4')
    fd = tempfile.TemporaryFile()
    print('test_1 =>',id(fd),id(fd.name))
    fd.write(b'test_1: my strings')
    #fd.seek(0)
    fd.flush()
    with open(fd.fileno()) as f:
        print('test_1: fd is opened and closed in with')
        print(f.readlines())


    fd = tempfile.TemporaryFile()
    print('test_2 =>',id(fd),id(fd.name))
    fd.write(b'test_2: my strings')

    # following will have:
    #
    #   OSError: [Errno 9] Bad file descriptor
    #
    # for any types of operations, e.gs;
    # fd.close(),  fd.flush()
    #
    # the name fd cannot be reused even in new instance
    #fd.seek(0)

    with open(fd.name) as f:
        print('test_2: fd is opened and closed in with')
        print(f.readlines())

import os
tmp1 = tempfile.TemporaryFile()
os.close(tmp1.fileno())                 # you may 


tmp2 = tempfile.TemporaryFile()
#print(id(tmp1), id(tmp2))
#del tmp1                       
print(id(tmp2))
tmp2.seek(0)

del tmp1
tmp2.write(b'good')
tmp2.seek(0)
