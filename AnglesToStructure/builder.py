"""
builds a molecule structure using its angles 

Boris Starkov(c)
"""


from math import cos, sin, sqrt

import Bio.PDB

from Bio.PDB import *


class Atom_Handy:
	'''
	Equivalent to pdb.Atom, but does special operations(moving & rotating) for building.
	'''

	def __init__(self, x, y, z, name):
		self.x = x
		self.y = y
		self.z = z
		self.name = name
	
	def move(self, dx, dy, dz):
		"""
		moves the atom along vector(dx, dy, dz)
		"""

		self.x += dx
		self.y += dy
		self.z += dz

	def rotate_over_arbitrary_axis(self, atom1, atom2, theta):
		"""
		rotates the atom itself about atom1 - atom2 line by angle theta.
		Positive angles are counter-clockwise looking down the axis toward the origin.
		The coordinate system is assumed to be right-hand.
		"""
		# let's make points from atoms
		p0 = [self.x, self.y, self.z]
		p1 = [atom1.x, atom1.y, atom1.z]
		p2 = [atom2.x, atom2.y, atom2.z]

		# Translate so axis is at origin
		p = [p0[i] - p1[i] for i in range(3)] 

		# Initialize point q
		q = [0.0,0.0,0.0]
		N = [p2[i] - p1[i] for i in range(3)]
		Nm = sqrt(N[0]**2 + N[1]**2 + N[2]**2)

		# Rotation axis unit vector
		n = [N[0] / Nm, N[1] / Nm, N[2] / Nm]

		# Matrix common factors     
		c = cos(theta)
		t = (1 - cos(theta))
		s = sin(theta)
		X = n[0]
		Y = n[1]
		Z = n[2]

		# Matrix 'M'
		d11 = t*X**2 + c
		d12 = t*X*Y - s*Z
		d13 = t*X*Z + s*Y
		d21 = t*X*Y + s*Z
		d22 = t*Y**2 + c
		d23 = t*Y*Z - s*X
		d31 = t*X*Z - s*Y
		d32 = t*Y*Z + s*X
		d33 = t*Z**2 + c

		#            |p[0]|
		# Matrix 'M'*|p[1]|
		#            |p[2]|
		q[0] = d11*p[0] + d12*p[1] + d13*p[2]
		q[1] = d21*p[0] + d22*p[1] + d23*p[2]
		q[2] = d31*p[0] + d32*p[1] + d33*p[2]

		# Translate axis and rotated point back to original location   

		result = [q[i] + p1[i] for i in range(3)]

		self.x = result[0]
		self.y = result[1]
		self.z = result[2]

	def to_atom(self):
		res = Atom.Atom(self.name, [self.x, self.y, self.z], 0, 0, '', self.name, None)
		return res


# calcs dihedral between for handy_atoms, not vectors

def calc_dihedral(a, b, c, d):
	v1 = a.to_atom().get_vector()
	v2 = b.to_atom().get_vector()
	v3 = c.to_atom().get_vector()
	v4 = d.to_atom().get_vector()
	return Bio.PDB.calc_dihedral(v1, v2, v3, v4)

# building the chain

def build_chain(angles):


	dist_n_ca = 1.46
	dist_ca_c = 1.51
	dist_c_n = 1.33

	mean_n_ca_c_angle = 1.941

	atoms = []

	atoms.append(Atom_Handy(0.0, 0.0, 0.0, 'N'))
	atoms.append(Atom_Handy(dist_n_ca, 0.0, 0.0, 'Ca'))
	atoms.append(Atom_Handy(dist_n_ca - cos(mean_n_ca_c_angle) * dist_ca_c, sin(mean_n_ca_c_angle) * dist_ca_c, 0.0, 'C'))

	for i in range(len(angles)):
		cur_angle = angles[i] # angle we want to make
		prev_atom = atoms[-1]
		if (i % 3 == 0):
			cur_atom = Atom_Handy(prev_atom.x + dist_c_n, prev_atom.y, prev_atom.z, 'N')
		elif (i % 3 == 1):
			cur_atom = Atom_Handy(prev_atom.x + dist_n_ca, prev_atom.y, prev_atom.z, 'Ca')
		else:
			cur_atom = Atom_Handy(prev_atom.x + dist_ca_c, prev_atom.y, prev_atom.z, 'C')

		cur_dihedral_angle = calc_dihedral(cur_atom, atoms[-1], atoms[-2], atoms[-3]) # current angle
		cur_atom.rotate_over_arbitrary_axis(atoms[-1], atoms[-2], -cur_angle + cur_dihedral_angle)
		atoms.append(cur_atom)
		
		# just checking everything is correct
		#rint(str(calc_dihedral(atoms[-1], atoms[-2], atoms[-3], atoms[-4])) + ' ' + str(cur_angle))
		#a1 = atoms[-1]
		#a2 = atoms[-2]
		#print(str(((a1.x - a2.x) ** 2 + (a1.y - a2.y) ** 2 + (a1.z - a2.z) ** 2) ** (1.0/2)) + ' ')
		#





build_chain([0.2, 0.566, 0.812, -1.7, 0.35])
























"""
################################################################################################################################
def move_atom(atom, ):

def rotate_atom(atom, ):
XYZ RotatePointAboutLine(XYZ p,double theta,XYZ p1,XYZ p2)

XYZ u,q1,q2;
double d;

/* Step 1 */
q1.x = p.x - p1.x;
q1.y = p.y - p1.y;
q1.z = p.z - p1.z;

u.x = p2.x - p1.x;
u.y = p2.y - p1.y;
u.z = p2.z - p1.z;
Normalise(&u);
d = sqrt(u.y*u.y + u.z*u.z);

/* Step 2 */
if (d != 0) {
  q2.x = q1.x;
  q2.y = q1.y * u.z / d - q1.z * u.y / d;
  q2.z = q1.y * u.y / d + q1.z * u.z / d;
} else {
  q2 = q1;
}

/* Step 3 */
q1.x = q2.x * d - q2.z * u.x;
q1.y = q2.y;
q1.z = q2.x * u.x + q2.z * d;

/* Step 4 */
q2.x = q1.x * cos(theta) - q1.y * sin(theta);
q2.y = q1.x * sin(theta) + q1.y * cos(theta);
q2.z = q1.z;

/* Inverse of step 3 */
q1.x =   q2.x * d + q2.z * u.x;
q1.y =   q2.y;
q1.z = - q2.x * u.x + q2.z * d;

/* Inverse of step 2 */
if (d != 0) {
  q2.x =   q1.x;
  q2.y =   q1.y * u.z / d + q1.z * u.y / d;
  q2.z = - q1.y * u.y / d + q1.z * u.z / d;
} else {
  q2 = q1;
}

/* Inverse of step 1 */
q1.x = q2.x + p1.x;
q1.y = q2.y + p1.y;
q1.z = q2.z + p1.z;
return(q1);
"""
