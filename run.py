from materials import material_properties
from flex_code import flex_code
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("num",help="./restart", type=str)
args=parser.parse_args()

material_properties=material_properties("start.xml")
run=flex_code(material_properties)

if args.num=="s":
    run.read_start("start.xml")
    run.write_to_restart("restart.xml")
else:
    run.read_restart("restart.xml")
run.flex_run()