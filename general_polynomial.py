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
				raise Exception('Exponent must be non-negative, exponent at index {idx} of exponents is negative.'.format(idx=idx))

		self.coef, self.exponents, self.variable_names = coef, exponents, variable_names

		if coef == 0:
			self.term_representation = '0'
		else:
			self.term_representation = '{coef}*'.format(coef=coef)+\
			'*'.join(
				[result for result in map(self.__variable_representation, variable_names, exponents) if result != '']
				)

		self.degree = max(exponents)

	@staticmethod
	def __variable_representation(variable_name: str, exponent: int):
		if exponent ==0:
			return ''
		elif exponent ==1:
			return variable_name
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
				continue
			elif value == 0:
				return 0
			else:
				total *= value**(self.exponents[idx])
		return total

	def __add__(self, other):
		if self.exponents == other.exponents:
			return GTerm(self.coef+other.coef, self.exponents, self.variable_names)
		else:
			#write a function that checks list of variables is identical and apply it to all operations
			#raise Exception('Terms must have the same list of exponents')
			return None

	def __sub__(self, other):
		if self.exponents == other.exponents:
			return GTerm(self.coef-other.coef, self.exponents, self.variable_names)
		else:
			#raise Exception('Terms must have the same list of exponents')
			return None

	def __mul__(self, other):
		new_exponents = tuple( [sum(val) for val in zip(self.exponents, other.exponents)] )
		return GTerm(self.coef*other.coef, new_exponents, self.variable_names)

	def __truediv__(self, other):
		new_exponents = tuple( [val[0]-val[1] for val in zip(self.exponents, other.exponents)] )
		for exponent in new_exponents:
			if exponent <0:
				raise Exception('Division by a lower degree term not allowed.')
		return GTerm(self.coef/other.coef, new_exponents, self.variable_names)

class GPolynomial:
	"""
	Build a polynomial from a list of term data [{'coeficient': float, 'exponents': tuple}]
	The only requirement is that the list of of exponents is fixed to the number of variables.
	"""
	def __init__(self, terms_data: list = [], variable_names: tuple = ()):
		self.terms_data = terms_data
		self.variable_names = variable_names 
		# TODO: implement inference for variables basis based on Terms data
		self.__build_terms()

	def __build_terms(self):
		self.terms=[]
		self.terms = [GTerm(term['coeficient'], term['exponents'], self.variable_names) for term in self.terms_data]

	def parse_polynomial(self, input_str:str):
		parser = GParser()
		terms_data, variable_names = parser.parse_input(input_str)
		return GPolynomial(terms_data, variable_names)

	def __repr__(self):

		base = '<GPolynomial :: '
		return  base + ' + '.join([str(term) for term in self.terms]) + '>'

	def __str__(self):

		return ' + '.join([str(term) for term in self.terms])

	def __call__(self, values: tuple) -> float:

		total = 0
		for term in self.terms:
			total += term(values)
		return total

	def __len__(self) -> int:
		# max of term degree
		return max([term.degree for term in self.terms])

	def reduce(self):
		"""reduces the polynomial by grouping and adding together similar terms. 
		It returns a new GPolynomial object with the reduced terms."""

		folded_terms_mapping = {}
		terms_to_skip = []
		for pointer1 in range(len(self.terms)):
			if pointer1 in terms_to_skip:
				continue

			temp_exponents = self.terms[pointer1].exponents
			temp_coef = self.terms[pointer1].coef
			folded_terms_mapping[pointer1] = {'terms_to_reduce': [], 'new_coef':temp_coef}

			for pointer2 in range(len(self.terms)):

				if pointer2 == pointer1:
					continue
				if temp_exponents == self.terms[pointer2].exponents:
					folded_terms_mapping[pointer1]['terms_to_reduce'].append(pointer2)
					folded_terms_mapping[pointer1]['new_coef'] += self.terms[pointer2].coef
					terms_to_skip.append(pointer2)

		new_terms_data = [{'coeficient': value['new_coef'], 'exponents': self.terms[key].exponents} for key, value in folded_terms_mapping.items() if value['new_coef'] !=0]

		return GPolynomial(new_terms_data, self.variable_names)

	def __add__(self, other):
		# new_variable_names = []
		# for i in self.variable_names + other.variable_names:
		# 	if i not in new_variable_names:
		# 		new_variable_names.append(i)
		if not self.variable_names == other.variable_names:
			raise Exception('Addition is not supported for GPolynomials with different variable name lists')

		new_polynomial = GPolynomial([], self.variable_names)
		new_polynomial.terms += self.terms + other.terms
		return new_polynomial.reduce()

	def __sub__(self, other):
		if not self.variable_names == other.variable_names:
			raise Exception('Substraction is not supported for GPolynomials with different variable name lists')

		for term in other.terms:
			term.coef = -term.coef
		return self + other

	def mult_term(self, other_term):
		if not self.variable_names == other_term.variable_names:
			raise Exception('Multiplive term must have the same variable list (basis?) as the polynomial')

		new_polynomial = GPolynomial([], self.variable_names)
		for term in self.terms:
			new_polynomial.terms.append(term*other_term)

		return new_polynomial.reduce()

	def __mul__(self, other):
		if not self.variable_names == other.variable_names:
			raise Exception('Multiplication is not supported for GPolynomials with different variable name lists')

		new_polynomial = GPolynomial([], self.variable_names)
		for term in other.terms:
			temp_p = self.mult_term(term)
			new_polynomial += temp_p

		return new_polynomial

	# TODO: implement division for GPolynomial

class GParser:
	def __init__(self):
		#input('write a polynomial expression')
		pass

	def parse_term(self, term_str: str):
		import re

		elements = term_str.split('*')
		variables_dict = {} #"{'<name>': exponent_value}"
		domain = []
		coef=1
		for el in elements:
			if re.match('-', el):
				coef *= -1
				el = el.replace('-', '')
			if re.match('[a-zA-Z]', el):
				if re.search('\^', el):
					var_name, exponent = el.split('^')[0], int(el.split('^')[1])
					if el.split('^')[0] not in domain:
						domain.append(el.split('^')[0])
				else:
					var_name, exponent = el, 1
					if el not in domain:
						domain.append(el)
				if variables_dict.get(var_name):
					variables_dict[var_name] += exponent
				else:
					variables_dict[var_name] = exponent

			else:
				if re.search('\^', el):
					coef *= float(el.split('^')[0])*el.split('^')[1]
				else:
					coef *= float(el)
		return {'coeficient': coef, 'variables_data': variables_dict, 'domain': domain}

	@staticmethod
	def __build_term_exponents(term_data_detailed, union_domain):
		exponents = [0]*len(union_domain)
		for idx, var in enumerate(union_domain):
			if term_data_detailed['variables_data'].get(var):
				exponents[idx] += term_data_detailed['variables_data'][var]
		return tuple(exponents)

	def parse_input(self, input_str: str):
		terms = input_str.split('+') # need to take into account minus sign, parentheses and dividion sign
		terms_pre_parse = []
		union_domain = []
		for term in terms:
			term_data_detailed = self.parse_term(term)
			union_domain.extend([var for var in term_data_detailed['domain'] if var not in  union_domain])
			terms_pre_parse.append(term_data_detailed)

		terms_data = []
		for term in terms_pre_parse:
			terms_data.append({'coeficient': term['coeficient'], 'exponents': self.__build_term_exponents(term, union_domain)})

		return terms_data, tuple(union_domain)

	def parse_raw_input(self, raw_str: str):
		# needs to parse parentheses ( ... ), division signe (...)/(...) etc. Each of those would be a polynomial to parse
		strings = input_str.split('/')


## temp tests
ex_parser = GParser()
input_str = '-5*x^2*y + 2*z'
print(ex_parser.parse_input(input_str))
poly_start = GPolynomial()
print(str(poly_start))
poly_parsed = poly_start.parse_polynomial("5*x^2*y+-1*y*z^3+6*x^2*y")
print(str(poly_parsed))
print(str(poly_parsed.reduce()))
print(poly_parsed.terms_data)
print(poly_parsed.variable_names)

print(GTerm(2,(1,2), ('x','y'))*GTerm(-3, (0,3), ('x','y')))
print(GTerm(2,(1,2), ('x','y'))/GTerm(-3, (2,1), ('x','y')))
		