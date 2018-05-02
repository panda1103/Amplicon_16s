#!/sas5/home/jiangdawei/python3/bin/python3
import argparse
import pandas as pd
import numpy as np
import matplotlib  
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import curve_fit
import os
import sys
import pdb
#def func(x, a, b):
#	return a * log(x) + b
#对数函数不合适

def func(x,a,b):
	return a*np.exp(b/x)

def col_map(arr):
	if len(arr) <= 8:
		col='Dark2'
	elif len(arr) <=10:
		col='tab10'
	elif len(arr) <=20:
		col='tab20'
	else:
		#col='gist_rainbow'
		#col='nipy_spectral'
		col='gist_ncar'
	return col

def plot_rare(rare, outfile):
	data = pd.read_table(rare)
	data=data.dropna()
	data=data.astype('float')
	data.index = data.pop("richness")
	#print(data)
	col = col_map(data.columns)
	for s in data.columns:
		y = data[s]
		#pdb.set_trace()
		#print(data.index)
		#print(y)
		popt, pcov = curve_fit(func, data.index, y)
		a = popt[0]
		b = popt[1]
		yvals=func(data.index, a, b)
		data[s] = yvals
	data.plot(colormap=col, xlim=(0,100))
	plt.subplots_adjust(left=0.1, right=0.80, top=0.9, bottom=0.2)
	plt.xlabel('Percent of reads')
	plt.ylabel('Number of OTU')
	#plt.legend(loc='upper right', fontsize=6)
	plt.legend(fontsize = 'xx-small', bbox_to_anchor=(1, 1))
	plt.title('Rarefaction curve')
	plt.savefig(outfile)

def plot_alpha(alpha, outfile):
	data = pd.read_table(alpha)
	data.index = data.pop("Sample")
	fig, axes = plt.subplots(3, 1)
	
	chao1 = data[['chao1']]
	#axes[0].set_title('chao1', fontsize=16)
	axes[0].set_ylabel('chao1', fontsize=10)
	chao1.plot.bar(subplots=False, color='r', ax=axes[0],sharex=True, sharey=False, fontsize=6, legend=False)

	shannon_e = data[['shannon_e']]
	#axes[1].set_title('shannon_e', fontsize=16)
	axes[1].set_ylabel('shannon_e', fontsize=10)
	shannon_e.plot.bar(subplots=False, color='g', ax=axes[1],sharex=True, sharey=False, fontsize=6, legend=False)

	simpson = data[['simpson']]
	#axes[2].set_title('simpson', fontsize=16)
	axes[2].set_ylabel('simpson', fontsize=10)
	#simpson.plot(figsize=(12, 6.5), colormap=col, subplots=False, ax=axes[2,0],sharex=True, sharey=False)
	simpson.plot.bar(subplots=False, color='b', ax=axes[2],sharex=True, sharey=False, fontsize=6, legend=False)
	plt.subplots_adjust(left=0.1, right=0.80, top=0.9, bottom=0.3)
	plt.savefig(outfile)

def plot_tax(tax, level, outfile):
	data = pd.read_table(tax)
	data.index = data.pop(level)
	data = data.T
	sam_num=data.shape[0]
	col = col_map(data.columns)
	if sam_num < 12:
		hight = sam_num
		length = sam_num
	elif sam_num <30:
		hight = 6
		length = sam_num/2
	elif sam_num <60:
		hight = 6
		length = sam_num/3
	else:
		hight = 6
		length = sam_num/4.5
	data.plot.bar(figsize=(length, hight), stacked=True, fontsize = 'small',  colormap=col)
	plt.legend(fontsize = 'xx-small', bbox_to_anchor=(1, 1))
	plt.subplots_adjust(left=0.1, right=0.80, top=0.9, bottom=0.35)
	#plt.show()
	plt.savefig(outfile)


def plot_beta_tree(tree, outfile, figpath):
	cmd = "java -Xmx4g -jar "+ figpath +" -graphic PDF -width 320 -height 320 " + tree + "  " + outfile
	os.system(cmd)

def main():
	parser = argparse.ArgumentParser(description='16s result visualization')
	parser.add_argument('--path', help='Analysis result path')
	basepath = sys.path[0]
	figpath = basepath + '/figtree.jar'
	args = parser.parse_args()
	inpath = args.path
	outpath = args.path
	#rare = inpath + "/alpha/rare_raw_part.txt"
	rare = inpath + "/alpha/rare_raw.txt"
	rareoutfile = outpath+ "/alpha/rarefaction_curve.pdf"
	plot_rare(rare, rareoutfile)

	alpha = inpath + "/alpha/alpha.txt"
	alphaoutfile = outpath+ "/alpha/alpha.pdf"
	plot_alpha(alpha, alphaoutfile)
	
	'''
	tree = inpath + "/beta/jaccard.tree"
	treeoutfile = outpath + "/beta/jaccard.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)

	tree = inpath + "/beta/jaccard_binary.tree"
	treeoutfile = outpath + "/beta/jaccard_binary.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)

	tree = inpath + "/beta/bray_curtis.tree"
	treeoutfile = outpath + "/beta/bray_curtis.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)

	tree = inpath + "/beta/bray_curtis_binary.tree"
	treeoutfile = outpath + "/beta/bray_curtis_binary.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)
	
	tree = inpath + "/beta/unifrac.tree"
	treeoutfile = outpath + "/beta/unifrac.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)

	tree = inpath + "/beta/unifrac_binary.tree"
	treeoutfile = outpath + "/beta/unifrac_binary.tree.pdf"
	plot_beta_tree(tree, treeoutfile, figpath)
	'''

	tax = inpath + "/tax/order_summary.txt"
	level = 'Order'
	taxoutfile = outpath + "/tax/order_summary.pdf"
	plot_tax(tax, level, taxoutfile)

	tax = inpath + "/tax/genus_summary.txt"
	level = 'Genus'
	taxoutfile = outpath + "/tax/genus_summary.pdf"
	plot_tax(tax, level, taxoutfile)

if __name__ == '__main__':
	main()
