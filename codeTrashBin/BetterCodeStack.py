import random

# Generate random IP Address
IP_address = ['.'.join((str(random.randint(1,254)) for t in range(4))) for s in range(10)]


