from __future__ import division
from copy import copy, deepcopy

class Term(object):
	def __init__(self, coef, order, var='x'):
		if order <0:
			raise Exception('Order must be non-negative')
		self.coef, self.order, self.var = coef, order, var

	def __repr__(self):
		base = 'Term({coef}, {order}) :: '.format(coef=self.coef, order=self.order)
		if self.coef==0:
			return base+'0'
		elif self.order==0:
			return base+'{coef}'.format(coef=self.coef)
		else:
			return base+'{coef}*{var}^{order}'.format(coef=self.coef, order=self.order, var=self.var)

	def __str__(self):
		if self.coef==0:
			return '0'
		elif self.order==0:
			return '{coef}'.format(coef=self.coef)
		else:
			return '{coef}*{var}^{order}'.format(coef=self.coef, order=self.order, var=self.var)

	def __call__(self, value):
		return self.coef*(value**self.order)

	def __add__(self, other):
		if self.order==other.order:
			return Term(self.coef+other.coef, self.order)
		else:
			raise Exception('Terms must be of same order')

	def __sub__(self, other):
		if self.order==other.order:
			return Term(self.coef-other.coef, self.order)
		else:
			raise Exception('Terms must be of same order')

	def __mul__(self, other):
		return Term(self.coef*other.coef, self.order + other.order)

	def __truediv__(self, other):
		if other.coef==0:
			raise Exception('Division by zero')
		if self.order >= other.order:
			return Term(self.coef/other.coef, self.order - other.order)
		else:
			raise Exception('numerator must have order greater or equal than denominator')

	def asPolynomial(self):
		trail = (self.order)*[0]
		trail.extend([self.coef])
		return Polynomial(trail)

class Polynomial(object):
	def __init__(self, coefs, roots=None, variable_name='x'):
		self.coefs = self.trim(coefs)
		self.roots = roots
		self.degree = len(coefs) -1
		self.var = variable_name
		self.update_terms()

	def __repr__(self):
		poly = ''
		for order in self.terms:
			if poly =='':
				poly = str(self.terms[self.degree-order])
			else:
				poly = poly + ' + '+str(self.terms[self.degree-order])
		return 'Polynomial({coefs}):: {algebra}'.format(coefs=self.coefs, algebra=poly)

	def __str__(self):
		poly = ''
		for order in self.terms:
			if poly =='':
				poly = str(self.terms[self.degree-order])
			else:
				poly = poly + ' + '+str(self.terms[self.degree-order])
		return '{algebra}'.format(coefs=self.coefs, algebra=poly)

	def __len__(self):
		return self.degree

	def __call__(self, value):
		res = 0
		for key in self.terms:
			res += self.terms[key](value)
		return res

	def update_terms(self):
		self.terms={}
		for order, coef in enumerate(self.coefs):
			self.terms[order]= Term(coef, order, var=self.var)

	def trim(self, coefs=None):
		trim = True
		if coefs:
			while trim:
				if coefs[len(coefs)-1]==0 and len(coefs)>1:
					del coefs[len(coefs)-1]
				else:
					trim = False
			return coefs
		else:
			while trim:
				if self.coefs[self.degree]==0 and len(self.coefs)>1:
					del self.coefs[self.degree]
					self.degree = len(self.coefs)-1
				else:
					trim = False
			self.update_terms()
			return self

	def __add__(self, other):
		delta = self.degree - other.degree
		if delta<0:
			self.coefs.extend(abs(delta)*[0])
		elif delta>0:
			other.coefs.extend(delta*[0])
		new_coefs = [x+y for x,y in zip(self.coefs, other.coefs)]
		new_coefs = self.trim(coefs=new_coefs)
		return Polynomial(new_coefs)

	def __sub__(self, other):
		delta = self.degree - other.degree
		if delta<0:
			self.coefs.extend(abs(delta)*[0])
		elif delta>0:
			other.coefs.extend(delta*[0])
		new_coefs = [x-y for x,y in zip(self.coefs, other.coefs)]
		new_coefs = self.trim(coefs=new_coefs)
		return Polynomial(new_coefs)

	def lead(self):
		return self.terms[self.degree]

	def mul_term(self, term):
		new_degree = self.degree + term.order
		new_terms = {}
		for order in self.terms:
			new_terms[order+term.order] = term*self.terms[order]
		new_coefs = (new_degree+1)*[0]
		for order in new_terms:
			new_coefs[order] = new_terms[order].coef
		return Polynomial(new_coefs)

	def __mul__(self, other):
		p = Polynomial([0])
		for order in other.terms:
			p_temp = self.mul_term(other.terms[order])
			p = p + p_temp.trim()
		p=p.trim()
		return p

	def __truediv__(self, other):
		quot = Polynomial([0])
		rem = copy(self)
		while (rem!=Polynomial([0])) and (len(rem) >= len(other)): 
			d = copy(other)
			t = rem.lead()/other.lead()
			quot = quot + t.asPolynomial()
			print('new quotient term: ', quot)
			rem = rem - d.mul_term(t)
			print('current reminder: ', rem)
		return (quot, rem)
	
	def derivative(self):
		derivative_coefs = []
		for order, coef in enumerate(self.coefs):
			derivative_coefs.append(order*coef)
		derivative_coefs = derivative_coefs[1:]
		return Polynomial(coefs=derivative_coefs, variable_name=self.var)

	@staticmethod
	def _next_guess(guess_k, guess_list, p, p_prime):
		a = p(guess_k)/p_prime(guess_k)
		b = sum([1/(guess_k - z) for z in guess_list if z != guess_k])
		w_k = a/(1 - a*b)
		return guess_k - w_k
	
	@staticmethod
	def _convergence_value(guess_list, p):
		return sum([abs(p(guess)) for guess in guess_list])
	
	@staticmethod
	def _round_complex(x, precision=4):
		return complex(round(x.real, precision),round(x.imag, precision))

	def find_roots(self, tolerance=10**-6, max_iters=100, precision=4):
		# https://en.wikipedia.org/wiki/Aberth_method
		from numpy.random import uniform
		p_prime = self.derivative()
		bound = 1 + max([coef/self.coefs[-1] for coef in self.coefs[:-1]]) # https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
		guess_list = uniform(-1, 1, self.degree) + 1.j * uniform(-1, 1, self.degree) # TODO: take into acocunt bound
		convergence_value = self._convergence_value(guess_list, self)
		iters = 0
		while (convergence_value > tolerance) & (iters < max_iters):
			iters+=1
			new_guess_list = []
			for guess in guess_list:
				new_guess_list.append(self._next_guess(guess, guess_list, self, p_prime))
			guess_list = new_guess_list
			convergence_value = self._convergence_value(guess_list, self)
		guess_list = [self._round_complex(x, precision) for x in guess_list]
		self.roots = guess_list
		print(f'Number of iterations was {iters}, while max iterations was {max_iters}. The covergence value is {convergence_value}'\
			.format(iters, max_iters, convergence_value))
		return guess_list






