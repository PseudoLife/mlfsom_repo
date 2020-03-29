from math import sin, cos, radians, sqrt

def CalcResolution(unit_cell,hkl):
	"""
	Calculates resolution of input hkl for space group
	unit_cell: tuple e.g. (77.25,77.25,38.66,90,90,90)
	hkl: tuple e.g. (3,2,2)
	"""
	a,b,c,alpha,beta,gamma = unit_cell
	alpha = radians(alpha)
	beta = radians(beta)
	gamma = radians(gamma)
	h,k,l = hkl

	V = a*b*c*sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*cos(alpha)*cos(beta)*cos(gamma))
	S11 = (b*c*sin(alpha))**2
	S22 = (a*c*sin(beta))**2
	S33 = (a*b*sin(gamma))**2
	S12 = a * b * c**2 * (cos(alpha)*cos(beta)-cos(gamma))
	S23 = a**2 * b * c * (cos(beta)*cos(gamma)-cos(alpha))
	S13 = a * b**2 * c * (cos(gamma)*cos(alpha)-cos(beta))

	inv_d_sqr = (1/V**2) * (S11*h**2 + S22*k**2 + S33*l**2 + 2*S12*h*k + 2*S23*k*l + 2*S13*h*l)
	d_space = inv_d_sqr**(-0.5)
	return d_space

"""
if __name__ == "__main__":
	unit_cell = (77.25,77.25,38.66,90,90,90)
	hkl = (3,2,2)
	unit_cell = (27.2800 ,31.9800 , 34.291, 88.530 , 108.500, 111.800)
	hkl = (3,2,2)
	print CalcResolution(unit_cell,hkl)
"""