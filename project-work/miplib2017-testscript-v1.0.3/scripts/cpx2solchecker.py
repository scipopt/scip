#!/usr/bin/env python
# encoding: utf-8
"""
cpx2solchecker.py

Convert a solution file from the CPLEX XML SOL format into the solution checker format
"""

from __future__ import print_function
import sys
import os
import gzip
from xml.etree.ElementTree import parse

usage="""
Convert a solution file from CPLEX XML SOL format into the solution checker format.

Usage:
{0} solution.sol

The mandatory argument 'solution.sol' is a solution in CPLEX XML SOL format.
"""


def cpx2sol(filename):
	"""
	Parse a Cplex XML solution file and prints the simplified SOL format

	=obj= objective value
	var1 value1
	var2 value2
	...

	Note: objective and variable values are never converted to and back from floats.
	They are kept as strings and moved from one format to the other.
	"""
	with open(filename) if not filename.endswith("gz") else gzip.open(filename, 'rb') as fp:
		doc = parse(fp).getroot()
		objVal = doc.find('.//header').get('objectiveValue')
		print("=obj=", objVal)
		for item in doc.findall('.//variables/variable'):
			print(item.get('name'), item.get('value'))


if __name__ == '__main__':
	if len(sys.argv) != 2:
		print(usage.format(sys.argv[0]))
		sys.exit()

	cpx2sol(sys.argv[1])
