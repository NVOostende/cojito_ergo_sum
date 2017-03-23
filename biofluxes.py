import numpy as np
import netCDF4 as nc 
import pyroms
import matplotlib.pylab as plt
from matplotlib.sankey import Sankey

class biofluxes():

	def __init__(self,pyroms_grid,file_dia_roms):

		self.grd = pyroms.grid.get_ROMS_grid(pyroms_grid)
		self.mask_subset = np.zeros(self.grd.hgrid.mask_rho.shape)

		self.file_dia_roms = file_dia_roms
		self.fact_units = 1035 * 86400 * 6.625 * 12# convert from mol N per kg/m per sec to mg C per m3/m per day

		self.dict_sankey = {'total_DIN':0,'Psm_DIN_uptake':0, 'Pmd_DIN_uptake':0, 'Plg_DIN_uptake':0,
		'Psm_exu_loss':0, 'Pmd_exu_loss':0,'Plg_exu_loss':0,
		'Psm_agg_loss':0, 'Pmd_agg_loss':0,'Plg_agg_loss':0,
		'Zsm_Psm_graz':0, 'Zmd_Pmd_graz':0, 'Zlg_Plg_graz':0,
		'Zsm_bact_graz':0, 'Zmd_Zsm_graz':0, 'Zlg_Zmd_graz':0,
		'Zsm_DOM_loss':0, 'Zmd_DOM_loss':0, 'Zlg_DOM_loss':0,
		'Zmd_det_loss':0, 'Zlg_det_loss':0, 'jhploss_Zmd':0, 'jhploss_Zlg':0,
		'Zdet_btf':0, 'det_btf': 0, 'NPP': 0, 'mesozooP': 0, 'hpp': 0,
		'pe-ratio': 0, 'z-ratio': 0, 'f-ratio': 0}

		self.list_nc_variables_dia_2d = ['npp_100',
		'prod_n_100_sm',
		'aggloss_n_100_sm',
		'zloss_n_100_sm',
		'prod_n_100_lg',
		'aggloss_n_100_lg',
		'zloss_n_100_lg',
		'prod_n_100_di',
		'aggloss_n_100_di',
		'zloss_n_100_di',
		'prod_n_100_smz',
		'ingest_n_100_smz',
		'zloss_n_100_smz',
		'hploss_n_100_smz',
		'prod_ndet_100_smz',
		'prod_n_100_mdz',
		'ingest_n_100_mdz',
		'zloss_n_100_mdz',
		'hploss_n_100_mdz',
		'prod_ndet_100_mdz',
		'prod_n_100_lgz',
		'ingest_n_100_lgz',
		'zloss_n_100_lgz',
		'hploss_n_100_lgz',
		'prod_ndet_100_lgz',
		'prod_n_100_bact',
		'zloss_n_100_bact',
		'uptake_din_100',
		'prod_mesozoo_100',
		'z_ratio_100',
		'pe_ratio_100',
		'f_ratio_100',
		'prod_don_100_smz',
		'prod_don_100_mdz',
		'prod_don_100_lgz',
		'prod_n_100_md',
		'aggloss_n_100_md',
		'zloss_n_100_md']
		

		return None

	def create_mask_analytic(self,imin,imax,jmin,jmax):
		self.mask_subset[jmin:jmax,imin:imax] = 1
		return None

	def compute_mean_of_2d_diag(self,varin):
		fid = nc.Dataset(self.file_dia_roms,'r')
		tmp = fid.variables[varin][:]
		out = (tmp * self.mask_subset * self.grd.hgrid.mask_rho * self.grd.hgrid.dx * self.grd.hgrid.dy).sum() / \
		(self.mask_subset * self.grd.hgrid.mask_rho * self.grd.hgrid.dx * self.grd.hgrid.dy).sum()
		fid.close()
		return out

	def populate_sankey_dict(self):
		# why is uptake_din_100 the same as npp_100? are they both in Nitrogen units?
		self.dict_sankey['total_DIN'] = self.compute_mean_of_2d_diag('uptake_din_100')
		self.dict_sankey['Psm_DIN_uptake'] = self.fact_units * self.compute_mean_of_2d_diag('prod_n_100_sm')
		self.dict_sankey['Pmd_DIN_uptake'] = self.fact_units * self.compute_mean_of_2d_diag('prod_n_100_md')
		self.dict_sankey['Plg_DIN_uptake'] = self.fact_units * self.compute_mean_of_2d_diag('prod_n_100_lg')
		self.dict_sankey['Psm_exu_loss'] = self.fact_units * 0.13 * self.compute_mean_of_2d_diag('prod_n_100_sm')
		self.dict_sankey['Pmd_exu_loss'] = self.fact_units * 0.13 * self.compute_mean_of_2d_diag('prod_n_100_md')
		self.dict_sankey['Plg_exu_loss'] = self.fact_units * 0.13 * self.compute_mean_of_2d_diag('prod_n_100_lg')
		self.dict_sankey['Psm_agg_loss'] = self.fact_units * self.compute_mean_of_2d_diag('aggloss_n_100_sm')
		self.dict_sankey['Pmd_agg_loss'] = self.fact_units * self.compute_mean_of_2d_diag('aggloss_n_100_md')
		self.dict_sankey['Plg_agg_loss'] = self.fact_units * self.compute_mean_of_2d_diag('aggloss_n_100_lg')
		self.dict_sankey['Zsm_Psm_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_sm')
		self.dict_sankey['Zmd_Pmd_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_md')
		self.dict_sankey['Zlg_Plg_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_lg')
		self.dict_sankey['Zsm_bact_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_bact')
		self.dict_sankey['Zmd_Zsm_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_smz')
		self.dict_sankey['Zlg_Zmd_graz'] = self.fact_units * self.compute_mean_of_2d_diag('zloss_n_100_mdz')
		self.dict_sankey['Zsm_DOM_loss'] = self.fact_units * self.compute_mean_of_2d_diag('prod_don_100_smz')
		self.dict_sankey['Zmd_DOM_loss'] = self.fact_units * self.compute_mean_of_2d_diag('prod_don_100_mdz')
		self.dict_sankey['Zlg_DOM_loss'] = self.fact_units * self.compute_mean_of_2d_diag('prod_don_100_lgz')
		self.dict_sankey['Zmd_det_loss'] = self.fact_units * self.compute_mean_of_2d_diag('prod_ndet_100_mdz')
		self.dict_sankey['Zlg_det_loss'] = self.fact_units * self.compute_mean_of_2d_diag('prod_ndet_100_lgz')
		self.dict_sankey['jhploss_Zmd'] = self.fact_units * self.compute_mean_of_2d_diag('hploss_n_100_mdz')
		self.dict_sankey['jhploss_Zlg'] = self.fact_units * self.compute_mean_of_2d_diag('hploss_n_100_lgz')
		self.dict_sankey['Zdet_btf'] = self.fact_units * (self.compute_mean_of_2d_diag('prod_ndet_100_mdz') + \
					       self.compute_mean_of_2d_diag('prod_ndet_100_lgz') )
		self.dict_sankey['det_btf'] = self.fact_units * self.compute_mean_of_2d_diag('ndet_btf')
		self.dict_sankey['NPP'] = self.compute_mean_of_2d_diag('npp_100')
		self.dict_sankey['mesozooP'] = self.compute_mean_of_2d_diag('prod_mesozoo_100')
		self.dict_sankey['hpp'] = self.fact_units * (self.compute_mean_of_2d_diag('hploss_n_100_smz') + \
		                          self.compute_mean_of_2d_diag('hploss_n_100_mdz') + self.compute_mean_of_2d_diag('hploss_n_100_lgz') )
		det_sink = self.dict_sankey['Zdet_btf'] + self.dict_sankey['Psm_agg_loss'] + self.dict_sankey['Pmd_agg_loss'] + \
			   self.dict_sankey['Plg_agg_loss']
		pe_ratio = np.max([0.,np.min([1., det_sink / self.compute_mean_of_2d_diag('npp_100')])])
		z_ratio = np.max([0.,np.min([1., self.compute_mean_of_2d_diag('mesozooprod_200') / self.compute_mean_of_2d_diag('npp_100')])])
		f_ratio = np.max([0.,np.min([1., self.compute_mean_of_2d_diag('uptake_no3_n2_100') / self.compute_mean_of_2d_diag('uptake_din_100')])])

		self.dict_sankey['pe-ratio'] = pe_ratio
		self.dict_sankey['z-ratio'] = z_ratio
		self.dict_sankey['f-ratio'] = f_ratio
		#...
		print self.dict_sankey
		return None

	def plot_sankey(self):
		fig = plt.figure(figsize=(15,10))
		ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], title='Biological C Flux Diagram')

		Nfluxes = self.dict_sankey
		sankey = Sankey(ax=ax, scale=1./Nfluxes['total_DIN'], offset=0.1, head_angle=140,
                format='%.0f', unit='')
		#prior 0
		sankey.add(patchlabel='', facecolor='lightgreen', rotation=180, trunklength=0.25,
		flows=[Nfluxes['total_DIN'],-Nfluxes['Psm_DIN_uptake'],-Nfluxes['Pmd_DIN_uptake'],-Nfluxes['Plg_DIN_uptake']],
		orientations=[0, -1, -1, -1],
		pathlengths=[0.25,0.5,0.26,0.52],
		labels=['DIN', None, None, None])
		#prior 1
		sankey.add(patchlabel='Psm', facecolor='lightgreen', rotation=180, trunklength=0.4,
		flows=[Nfluxes['Psm_DIN_uptake'],-Nfluxes['Psm_exu_loss'],-Nfluxes['Psm_agg_loss'],-Nfluxes['Zsm_Psm_graz']],
		orientations=[1, 1, 1, 0],
		labels=['Psm\njprod', 'exu', 'agg', None],
		pathlengths=[0.25,0.1,0.7,0.5],
		prior=0, connect=(1,0))
		#prior 2
		sankey.add(patchlabel='Zsm', facecolor='lightblue', trunklength=0.4,
		flows=[Nfluxes['Zsm_Psm_graz'],Nfluxes['Zsm_bact_graz'],-Nfluxes['Zsm_DOM_loss'],-Nfluxes['Zmd_Zsm_graz']],
		orientations=[0, 1, 1, -1],
		labels=['Zsm_jprod', 'bact', 'DOC', ''],
		pathlengths=[0.1,0.5,0.1,0.12],
		prior=1, connect=(3,0))
		#prior 3
		sankey.add(patchlabel='Pmd', facecolor='lightgreen', rotation=180, trunklength=1.0,
		flows=[Nfluxes['Pmd_DIN_uptake'],-Nfluxes['Pmd_exu_loss'],-Nfluxes['Pmd_agg_loss'],-Nfluxes['Zmd_Pmd_graz']],
		orientations=[1, 1, 1, 0],
		labels=['Pmd\njprod', 'exu', 'agg', None],
		pathlengths=[0.55,0.1,1.35,0.50],
		prior=0, connect=(2,0))
		#prior 4
		sankey.add(patchlabel='Zmd', facecolor='lightblue', trunklength=.6,
		flows=[Nfluxes['Zmd_Pmd_graz'],Nfluxes['Zmd_Zsm_graz'],-Nfluxes['Zmd_DOM_loss'],-Nfluxes['Zmd_det_loss'],-Nfluxes['Zlg_Zmd_graz'],-Nfluxes['jhploss_Zmd']],
		orientations=[0, 1, 1, 1, -1, -1],
		labels=['Zmd_jprod', '', 'DOC', 'det', '','hp'],
		pathlengths=[0.4,0.1,0.1,1.0,0.1,1.0],
		prior=2, connect=(3,1))
		#prior 5
		sankey.add(patchlabel='Plg', facecolor='lightgreen', rotation=180, trunklength=1.85,
		flows=[Nfluxes['Plg_DIN_uptake'],-Nfluxes['Plg_exu_loss'],-Nfluxes['Plg_agg_loss'],-Nfluxes['Zlg_Plg_graz']],
		orientations=[1, 1, 1, 0],
		labels=['Plg\njprod', 'exu', 'agg', None],
		pathlengths=[0.55,0.1,2.2,0.52],
		prior=0, connect=(3,0))
		#prior 6
		sankey.add(patchlabel='Zlg', facecolor='lightblue', trunklength=1.,
		flows=[Nfluxes['Zlg_Plg_graz'],Nfluxes['Zlg_Zmd_graz'],-Nfluxes['Zlg_DOM_loss'],-Nfluxes['Zlg_det_loss'],-Nfluxes['jhploss_Zlg']],
		orientations=[0, 1, 1, 1, -1],
		labels=['Zmd_jprod', '', 'DOC', 'det', 'hp'],
		pathlengths=[1.25,0.1,0.1,1.8,0.5],
		prior=4, connect=(4,1))
 
		sankey.finish()
 
		fig.text(0.60,0.75, 'NPP: {:.0f}'.format(Nfluxes['NPP'])+' mg C $\mathregular{m^{-2}}$ $\mathregular{d^{-1}}$',\
			 fontsize=12, fontweight='bold')
		fig.text(0.22,0.45, 'mesozooP:\n {:.0f}'.format(Nfluxes['mesozooP'])+' mg C $\mathregular{m^{-2}}$ $\mathregular{d^{-1}}$',\
			 fontsize=12, fontweight='bold')
		fig.text(0.30,0.25, 'Benthic flux: {:.0f}'.format(Nfluxes['det_btf'])+' mg C $\mathregular{m^{-2}}$ $\mathregular{d^{-1}}$',\
			 fontsize=12, fontweight='bold')
		fig.text(0.30,0.85, 'Pelagic flux: {:.0f}'.format(Nfluxes['hpp'])+' mg C $\mathregular{m^{-2}}$ $\mathregular{d^{-1}}$',\
			 fontsize=12, fontweight='bold')
		fig.text(0.80,0.87, 'pe-ratio: {:4.2f}'.format(Nfluxes['pe-ratio']), fontsize=12) #need to be checked, not a 0-1 ratio
		fig.text(0.80,0.84, 'z-ratio: {:6.2f}'.format(Nfluxes['z-ratio']), fontsize=12)   #need to be checked
		fig.text(0.80,0.81, 'f-ratio: {:6.2f}'.format(Nfluxes['f-ratio']), fontsize=12)   #need to be checked
		plt.show()
		return None






#-------------- test -----------------------

#bf1 = biofluxes('CCS1','/Volumes/P4/workdir/raphael/test_cojito/CCS1-RD.NVOcobalt31S_dia_00053.nc')
bf1 = biofluxes('CCS1','/Users/nicolasvanoostende/Dropbox/CCS_ROMS/new_online_diags/CCS1-RD.NVOcobalt31S_dia_00053.nc')
bf1.create_mask_analytic(130,180,300,450) # coastal north CCS
#bf1.create_mask_analytic(10,60,230,290) # offshore central CCS
#bf1.create_mask_analytic(100,123,230,290) # coastal central CCS
bf1.populate_sankey_dict()
bf1.plot_sankey()
