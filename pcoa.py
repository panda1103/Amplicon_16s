#!/sas5/home/jiangdawei/python3/bin/python3
# ----------------------------------------------------------------------------
# The PCoA method is modified from scikit-bio.
# jiangdawei@icarbonx.com
# 2018-1-31
# ----------------------------------------------------------------------------
import argparse
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from functools import partial
import warnings

class PCoA():
	r"""Perform Principal Coordinate Analysis.

	   If the distance is not euclidean (for example if it is a
	   semimetric and the triangle inequality doesn't hold),
	   negative eigenvalues can appear. There are different ways
	   to deal with that problem (see Legendre & Legendre 1998, \S
	   9.2.3), but none are currently implemented here.

	   However, a warning is raised whenever negative eigenvalues
	   appear, allowing the user to decide if they can be safely
	   ignored.
	"""
	short_method_name = 'PCoA'
	long_method_name = 'Principal Coordinate Analysis'

	def __init__(self, matrix, ids):
		self.dm = np.asarray(matrix, dtype=np.float64)
		self.ids = ids
		self._pcoa()

	def _pcoa(self):
		E_matrix = self._E_matrix(self.dm)

		# If the used distance was euclidean, pairwise distances
		# needn't be computed from the data table Y because F_matrix =
		# Y.dot(Y.T) (if Y has been centred).
		F_matrix = self._F_matrix(E_matrix)

		# If the eigendecomposition ever became a bottleneck, it could
		# be replaced with an iterative version that computes the
		# largest k eigenvectors.
		eigvals, eigvecs = np.linalg.eigh(F_matrix)

		# eigvals might not be ordered, so we order them (at least one
		# is zero). cogent makes eigenvalues positive by taking the
		# abs value, but that doesn't seem to be an approach accepted
		# by L&L to deal with negative eigenvalues. We raise a warning
		# in that case. First, we make values close to 0 equal to 0.
		negative_close_to_zero = np.isclose(eigvals, 0)
		#print(eigvals)
		eigvals[negative_close_to_zero] = 0
		#print(eigvals)
		if np.any(eigvals < 0):
		    warnings.warn(
		        "The result contains negative eigenvalues."
		        " Please compare their magnitude with the magnitude of some"
		        " of the largest positive eigenvalues. If the negative ones"
		        " are smaller, it's probably safe to ignore them, but if they"
		        " are large in magnitude, the results won't be useful. See the"
		        " Notes section for more details. The smallest eigenvalue is"
		        " {0} and the largest is {1}.".format(eigvals.min(),
		                                              eigvals.max()),
		        RuntimeWarning
		        )
		idxs_descending = eigvals.argsort()[::-1]
		self.eigvals = eigvals[idxs_descending]
		self.eigvecs = eigvecs[:, idxs_descending]

	def scores(self):
		"""Compute coordinates in transformed space.
		"""
		# Scale eigenvalues to have lenght = sqrt(eigenvalue). This
		# works because np.linalg.eigh returns normalized
		# eigenvectors. Each row contains the coordinates of the
		# objects in the space of principal coordinates. Note that at
		# least one eigenvalue is zero because only n-1 axes are
		# needed to represent n points in an euclidean space.

		# If we return only the coordinates that make sense (i.e., that have a
		# corresponding positive eigenvalue), then Jackknifed Beta Diversity
		# won't work as it expects all the OrdinationResults to have the same
		# number of coordinates. In order to solve this issue, we return the
		# coordinates that have a negative eigenvalue as 0
		num_positive = (self.eigvals >= 0).sum()
		eigvecs = self.eigvecs
		eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
		eigvals = self.eigvals
		eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)

		coordinates = eigvecs * np.sqrt(eigvals)

		proportion_explained = eigvals / eigvals.sum()
		
		return (eigvals, coordinates, proportion_explained, self.ids)
		r'''return OrdinationResults(eigvals=eigvals, site=coordinates,
		                         proportion_explained=proportion_explained,
		                         site_ids=self.ids)
		'''

	@staticmethod
	def _E_matrix(distance_matrix):
		"""Compute E matrix from a distance matrix.

		Squares and divides by -2 the input elementwise. Eq. 9.20 in
		Legendre & Legendre 1998."""
		return distance_matrix * distance_matrix / -2

	@staticmethod
	def _F_matrix(E_matrix):
		"""Compute F matrix from E matrix.

		Centring step: for each element, the mean of the corresponding
		row and column are substracted, and the mean of the whole
		matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
		row_means = E_matrix.mean(axis=1, keepdims=True)
		col_means = E_matrix.mean(axis=0, keepdims=True)
		matrix_mean = E_matrix.mean()
		return E_matrix - row_means - col_means + matrix_mean

def _validate_plot_axes(coord_matrix, axes):
	"""Validate `axes` against coordinates matrix."""
	num_dims = coord_matrix.shape[0]
	if num_dims < 3:
		raise ValueError("At least three dimensions are required to plot "
						"ordination results. There are only %d "
						"dimension(s)." % num_dims)
		if len(axes) != 3:
			raise ValueError("axes must contain exactly three elements (found "
		                     "%d elements)." % len(axes))
		if len(set(axes)) != 3:
			raise ValueError("The values provided for axes must be unique.")

		for idx, axis in enumerate(axes):
		    if axis < 0 or axis >= num_dims:
		        raise ValueError("axes[%d] must be >= 0 and < %d." %
									(idx, num_dims))

def _get_plot_point_colors(df, column, ids, cmap):
	if((df is None and column is not None) or (df is not None and
		column is None)):
		raise ValueError("Both df and column must be provided, or both "
				"must be None.")
	elif df is None and column is None:
		point_colors = None
	else:
		if column not in df:
			raise ValueError("Column '%s' not in data frame." % column)
		col_vals = df.loc[ids, column]
		if col_vals.isnull().any():
			raise ValueError("One or more IDs in the ordination results "
								"are not in the data frame, or there is "
								"missing data in the data frame's '%s' "
								"column." % column)
		category_to_color = None
		try:
			point_colors = col_vals.astype(float)
		except ValueError:
			categories = col_vals.unique()
			cmap = plt.get_cmap(cmap)
			category_colors = cmap(np.linspace(0, 1, len(categories)))
			category_to_color = dict(zip(categories, category_colors))
			point_colors = col_vals.apply(lambda x: category_to_color[x])
		point_colors = point_colors.tolist()
		#print(point_colors)
	return point_colors

def plot(coordinates, ids, outfile, df=None, column=None, axes=(0, 1, 2), axis_labels=('PC1', 'PC2', 'PC3'),
		title='Principal Coordinate Analysis', cmap='RdYlBu', s=30):
	coord_matrix = coordinates.T
	_validate_plot_axes(coord_matrix, axes)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	point_colors = _get_plot_point_colors(df, column, ids, cmap)
	xs = coord_matrix[axes[0]]
	ys = coord_matrix[axes[1]]
	zs = coord_matrix[axes[2]]
	scatter_fn = partial(ax.scatter, xs, ys, zs, s=s)
	ax.scatter(xs, ys, zs, s=s, c=point_colors)
	if point_colors is None:
		point_colors = 'b'
	plot = scatter_fn(c=point_colors, cmap=cmap)
	for x, y, z, l in zip(xs, ys, zs, ids):
		#ax.text(x, y, z, l, fontsize='xx-small')
		ax.text(x, y, z, l, fontsize=3.5)
	ax.set_xlabel(axis_labels[0])
	ax.set_ylabel(axis_labels[1])
	ax.set_zlabel(axis_labels[2])
	ax.set_title(title)
	plt.savefig(outfile)
	#plt.show()

def main():

	parser = argparse.ArgumentParser(description='Principal Coordinate Analysis(PCoA)')
	parser.add_argument('--distance_file')
	parser.add_argument('--distance_type', default='unifrac')
	parser.add_argument('--group_file')
	parser.add_argument('--group_lable', default=None)
	parser.add_argument('--outfile_prefix', default='pcoa')
	args = parser.parse_args()
	distance = args.distance_file
	distance_type = args.distance_type

	if args.group_file:
		group_file = args.group_file
		df = pd.read_table(group_file)
		df.index = df.pop('#sample')
	if args.group_lable:
		column = args.group_lable
	
	outfile_prefix = args.outfile_prefix
	outfile_result = outfile_prefix + '.pcoa.csv'
	outfile_result2 = outfile_prefix + '.pcoa.proportion_explained.csv'
	outfile_pdf = outfile_prefix + '.pcoa.pdf'
	
	data = pd.read_table(distance)
	data.index = data.pop(distance_type)
	matrix = np.array(data)
	ids = data.index
	pcoa_result = PCoA(matrix, ids)
	eigvals, coordinates, proportion_explained, ids = pcoa_result.scores()
#	print(test.scores())
	pcoa = pd.DataFrame(coordinates,index=ids)
	pcoa = pcoa.iloc[0:,0:3]
	pcoa.rename(columns={0:'PC1', 1:'PC2', 2:'PC3'}, inplace=True)
	pcoa.to_csv(outfile_result)
	pcoa_labels=('PC1', 'PC2', 'PC3')
	proportion_explained = map(lambda x:'%.2f' % x, [100*y for y in proportion_explained])
	proportion_explained = list(proportion_explained)[0:3]
	axis_labels=[]
	for x,y in zip(pcoa_labels, proportion_explained):
		axis_labels.append(x+'('+y+'%)')
	#print(axis_labels)
	plot(coordinates=coordinates, ids=ids, outfile=outfile_pdf, axis_labels=axis_labels, df=df, column=column)
	proportion_explained = pd.DataFrame(proportion_explained, index=pcoa_labels, columns=['Proportion explained(%)'])
	proportion_explained.to_csv(outfile_result2)

if __name__=='__main__':
	main()
