from __future__ import division
from copy import copy, deepcopy
import typing

class GTerm:
	""" 
	Generalized Polynomial term. Defined by a coefiecient and a list of exponents. 
	example 1: {'coeficient': 4, 'exponents': [1,3,2]} represents 4*x*y^3*z^2
	example 2: {'coeficient': 5, 'exponents': [1,0,2]} represents 5*x*z^2
	"""
	def __init__(self, coef: float, exponents: tuple, variable_names: tuple):

		for idx, exponent in enumerate(exponents):
			if exponent <0:
				raise Exception('Exponent must be non-negative, exponent at index {idx} is negative.'.format(idx=idx))

		self.coef, self.exponents, self.variable_names = coef, exponents, variable_names

		self.term_representation = '{coef}*'.format(coef=coef)+\
		'*'.join(
			[result for result in map(self.__variable_representation, variable_names, exponents) if result != '']
			)
		self.degree = sum(exponents)

	@staticmethod
	def __variable_representation(variable_name: str, exponent: int):
		if exponent ==0:
			return ''
		else:
			return '{var}^{exponent}'.format(var=variable_name, exponent=exponent)

	
	def __repr__(self):
		base = '<GTerm :: '
		return base + self.term_representation + '>'

	def __str__(self):
		return self.term_representation

	def __call__(self, values: tuple) -> float:

		#raise if values is not same as number of variables

		total = self.coef
		for idx, value in enumerate(values):
			if self.exponents[idx] == 0:
				total *=1
			elif value == 0:
				return 0
			else:
				total *= value**(self.exponents[idx])
		return total

	@staticmethod
	def __absolute_diff(a, b):
		return abs(a-b)

	def __add__(self, other):
		if sum(map(self.__absolute_diff, self.exponents, other.exponents)) ==0:
			return GTerm(self.coef+other.coef, self.exponents, self.variable_names)
		else:
			#write a function that checks list of variables is identical and apply it to all operations
			raise Exception('Terms must have the same list of exponents')

	def __sub__(self, other):
		if sum(map(self.__absolute_diff, self.exponents, other.exponents)) ==0:
			return GTerm(self.coef-other.coef, self.exponents, self.variable_names)
		else:
			raise Exception('Terms must have the same list of exponents')

	@staticmethod
	def __sum_exponents(a, b):
		return a+b

	def __mul__(self, other):

		new_exponents = tuple( map(self.__sum_exponents, self.exponents, other.exponents) )

		return GTerm(self.coef*other.coef, new_exponents, self.variable_names)

class GPolynomial:
	"""
	Build a polynomial from a list of arebitrary term data [{'coeficient': float, 'exponents': tuple}]
	The only requirement is that the list of of exponents is fixed to the number of variables.
	"""
	def __init__(self, terms_data: list, variable_names: tuple):
		self.terms_data = terms_data
		self.variable_names = variable_names
		self.__build_terms()

	def __build_terms(self):
		self.terms=[]
		self.terms = [GTerm(term['coeficient'], term['exponents'], self.variable_names) for term in self.terms_data]

	def __repr__(self):

		base = '<GPolynomial :: '
		
		return  base + ' + '.join([str(term) for term in self.terms]) + '>'

	def __str__(self):

		return ' + '.join([str(term) for term in self.terms])

	def __call__(self, values):

		total = 0
		for term in self.terms:
			total += term(values)
		return total

	def __len__(self):
		# max sum of term exponent (i.e., term degree)
		return max([term.degree for term in self.terms])






## temp tests
term = GTerm(4, (2,0,3), ('x','y','z'))
print(repr(term))
print(term)
print(term((-3,10,1)))
term2 = GTerm(3, (1,1,0), ('x','y','z'))
print(str(term2))

print(term*term2)
#print(term + term2)
poly = GPolynomial([{'coeficient':5, 'exponents': (5,2)}, {'coeficient':-3, 'exponents': (0,2)}], ('x', 'y'))
print(repr(poly))
print(str(poly))
print(poly((1,5)))
print(len(poly))
		