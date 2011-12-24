'''
Created on Dec 24, 2011

@author: cardona
'''

import functools

class memoize(object):
	"""Decorator that caches a function's return value each time it is called.
	If called later with the same arguments, the cached value is returned, and
	not re-evaluated.
	The cache is stored as a property of the instance decorated, so that it is
	very easy to clear it.
	"""
	def __init__(self, func):
		print "Init"
		self.func = func
		#self.cache = {}
		
	def __call__(self, *args):
		print "Call"
		if not self.func in self.cache:
			 self.cache[self.func] = {}
		try:
			return self.cache[self.func][args]
		except KeyError:
			value = self.func(*args)
			self.cache[self.func][args] = value
			return value
		except TypeError:
			# uncachable -- for instance, passing a list as an argument.
			# Better to not cache than to blow up entirely.
			return self.func(*args)
	  
	def __repr__(self):
		"""Return the function's docstring."""
		return self.func.__doc__
  
	def __get__(self, obj, objtype):
		"""Support instance methods."""
		print "Get"
		self.cache = obj.cache
		fn = functools.partial(self.__call__, obj)
		return fn
  
	def _echo(self,*args):
		 print self.cache
		 print args
