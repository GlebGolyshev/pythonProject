import argparse
parser = argparse.ArgumentParser()
parser.add_argument("srcf", help="input file name to check")
parser.add_argument("srct", help="input unified table location")
parser.add_argument("Dir", help="output directory")
args = parser.parse_args()

a = args.srcf
b = args.srct
c = args.Dir
print(a, 'file', b, 'folder', c, 'dir')
# parser = argparse.ArgumentParser()
# parser.add_argument("square", type=int,
#                     help="display a square of a given number")
# parser.add_argument("-v", "--verbosity", type=int,
#                     help="increase output verbosity")
# args = parser.parse_args()
# answer = args.square**2
# if args.verbosity == 2:
#     print(f"the square of {args.square} equals {answer}")
# elif args.verbosity == 1:
#     print(f"{args.square}^2 == {answer}")
# else:
#     print(answer)