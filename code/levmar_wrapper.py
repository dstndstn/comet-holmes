import levmar


class Optimizer(object):

	# func, jacobian, init, measurements, extra,
	# low, high, offset, scale
	def __init__(self, func, jacobian=None, initial=None,
				 measurements=None, extra=None, low=None,
				 high=None, offset=None, scale=None):
		self.func = func
		self.jacobian = jacobian
		self.initial = initial
		self.measurements = measurements
		self.extra = extra
		self.low = low
		self.high = high
		self.offset = offset
		self.scale = scale
		# Step size of approximate Jacobian.
		self.jacobian_delta = None

	def set_initial(self, init):
		self.initial = init

	def set_low(self, low):
		self.low = low

	def set_high(self, high):
		self.high = high

	def set_offset(self, offset):
		self.offset = offset

	def scaleoffset(self, params):
		# scale params
		params = params[:]
		if self.offset:
			params = [p - o for (p,o) in zip(params, self.offset)]
		if self.scale:
			params = [p / s for (p,s) in zip(params, self.scale)]
		return params

	def invscaleoffset(self, params):
		# scale params
		params = params[:]
		if self.scale:
			params = [p * s for (p,s) in zip(params, self.scale)]
		if self.offset:
			params = [p + o for (p,o) in zip(params, self.offset)]
		return params

	def objective(self, params, extra):
		#print 'before scale/offset:', params
		params = self.invscaleoffset(params)
		#print 'after scale/offset:', params
		model = self.func(params, self.measurements, extra)
		return model

	@staticmethod
	def theobjective(params, target, extra):
		(me, extra) = extra
		return me.objective(params, extra)

	def optimize(self):
		extra = (self, self.extra)

		lower = self.scaleoffset(self.low)
		upper = self.scaleoffset(self.high)

		iterations = 10

		opts = list(levmar.DEFAULT_OPTS)
		if not self.jacobian:
			if self.jacobian_delta:
				d = self.jacobian_delta
			else:
				d = levmar.DIFF_DELTA
			opts.append(d)

		# ||J'e||_inf
		opts[1] = 1e-10
		# Dp, change in parameter values.
		opts[2] = 1e-6
		# ||e||_2, norm of error.
		# opts[3] = 1e-3

		if self.initial is None:
			init = (0,) * len(self.scale)
		else:
			init = self.scaleoffset(self.initial)

		#print 'initial point:', init
		#print 'low:', lower
		#print 'high:', upper
		#print 'scale:', self.scale
		#print 'offset:', self.offset

		(result, iterations, run_info) = levmar.ddif_bc(
			Optimizer.theobjective,
			init,
			self.measurements,
			lower,
			upper,
			iterations,
			opts=opts,
			data=extra)

		if result:
			result = self.invscaleoffset(result)
		return result


