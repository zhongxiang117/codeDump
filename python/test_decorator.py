import time

def func_timing(func):

    """
    Outputs the time a function takes to execute.
    """

    def wrapper():
        t1 = time.time()
        func()
        t2 = time.time()
        return "Executing time: " + str((t2 - t1)) + "\n"
    return wrapper


@func_timing
def my_func():
    nmlist = [i for i in range(1000)]
    print("\nSum of all the numbers: " + str((sum(nmlist))))


print(my_func())




