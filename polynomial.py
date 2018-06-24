from __future__ import division
from copy import copy, deepcopy

class Term(object):
	def __init__(self, coef, order):
		if order <0:
			raise Exception('Order must be non-negative')
		self.coef, self.order = coef, order

	def __repr__(self):
		base = 'Term({coef}, {order})'.format(coef=self.coef, order=self.order)
		if self.coef==0:
			return base+':: 0'
		elif self.order==0:
			return base+':: {coef}'.format(coef=self.coef)
		else:
			return base+':: {coef}*x^{order}'.format(coef=self.coef, order=self.order)

	def __str__(self):
		if self.coef==0:
			return '0'
		elif self.order==0:
			return '{coef}'.format(coef=self.coef)
		else:
			return '{coef}*x^{order}'.format(coef=self.coef, order=self.order)

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
	def __init__(self, coefs):
		self.coefs = self.trim(coefs)
		self.degree = len(coefs) -1
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
			self.terms[order]= Term(coef, order)

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
			print 'new quotient term: ', quot
			rem = rem - d.mul_term(t)
			print 'current reminder: ', rem
		return (quot, rem)

