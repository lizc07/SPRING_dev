#!/home/Public/BioSoft/anaconda2/bin/python
# lizc@20180822
from argparse import ArgumentParser
from spring_helper import *
from doublet_detector import *
from collections import defaultdict
import csv
def r_equal(a, b):
  if len(a) == len(b) and all([a[i] == b[i] for i in range(0,len(a))]):
    return True
  else:
    return False

def r_match(a, b):
  [b.index(x) if x in b else None for x in a]

def subset_data_from_csv(csvdata, rownames = None, colnames = None): # data is from read_csv
  if not rownames:
    rownames = csvdata['rownames']
  if not colnames:
    colnames = csvdata['colnames']
    
  row_ind = [csvdata['rownames'].tolist().index(i) for i in rownames]
  col_ind = [csvdata['colnames'].tolist().index(i) for i in colnames]
  retdata = {
    'rownames': np.array(rownames),
    'colnames': np.array(colnames),
    'data'    : np.array(csvdata['data'][row_ind,:][:,col_ind])
    }
  return retdata
  
def read_csv(file, is_numeric = False):
  if not file:
    return {}
  count = 0
  rownames = []
  colnames = []
  data = []
  for row in csv.reader(open(file,'r')):
    if count == 0:
      colnames = row[1:]
    else:
      rownames.append(row[0])
      if is_numeric:
        data.append([float(i) for i in row[1:]])
      else:
        data.append([i for i in row[1:]])  
    count += 1
  return {'rownames' : np.array(rownames), 'colnames' : np.array(colnames), 'data' : np.array(data)}
# tsne_data = read_csv('Fetal_liver_10x_tsne_coord.csv', True)
# meta_data = read_csv('Fetal_liver_10x_cluster_info.csv', False)

def generate_data(input_path, output_path, meta_data, tsne_data):
  print 'loading UMI matrix ...'
  if os.path.isfile(input_path + '/matrix.npz'):
    E = scipy.sparse.load_npz(input_path + '/matrix.npz')
  else:
    E = scipy.io.mmread(input_path + '/matrix.mtx').T.tocsc()
    #scipy.sparse.save_npz(input_path + '/matrix.npz', E, compressed=True)
  print 'UMI matrix is loaded!'
  print E.shape
  
  print 'loading gene list ...'
  gene_list = np.array(load_genes(input_path + '/genes.tsv', delimiter='\t', column=1))
  print 'loading barcode list ...'
  barcode_list = np.array([bar.strip('\n') for bar in open(input_path + '/barcodes.tsv')])
  # sometimes, the barcode has '-1' suffix, so we will keep first 16 characters only for 10X data.
  barcode_list = [b[0:16] for b in barcode_list]
  cell_barcode = barcode_list
  cell_groupings = {}

  # load meta.data and tsne.data
  if meta_data:
    print 'loading meta data ...'
    meta_data = read_csv(meta_data, False)
    cell_barcode = np.intersect1d(cell_barcode, meta_data['rownames'])
  if tsne_data:
    print 'loading tsne data ...'
    tsne_data = read_csv(tsne_data, True)
    cell_barcode = np.intersect1d(cell_barcode, tsne_data['rownames'])
  
  # subset datasets with cell_barcode retained for subsequent analysis
  if not r_equal(cell_barcode, barcode_list):
    print 'subsetting UMI matrix...'
    E = E[[barcode_list.index(i) for i in cell_barcode],:]
    if meta_data:
      print 'subsetting meta data ...'
      meta_data = subset_data_from_csv(meta_data, cell_barcode.tolist())
      for i in range(0,len(meta_data['colnames'])):
        cell_groupings[meta_data['colnames'][i]] = meta_data['data'].T[i].tolist()
    if tsne_data:
      print 'subsetting tsne data ...'
      tsne_data = subset_data_from_csv(tsne_data, cell_barcode.tolist())
    
  total_counts = E.sum(1).A.squeeze()
  E = tot_counts_norm(E)[0]

  main_spring_dir = output_path+'/'
  if not os.path.exists(main_spring_dir):
    print 'making output_path: ' + main_spring_dir
    os.makedirs(main_spring_dir)
  
  print 'saving genes.txt and total_counts.txt ...'
  np.savetxt(main_spring_dir + 'genes.txt', gene_list, fmt='%s')
  np.savetxt(main_spring_dir + 'total_counts.txt', total_counts)
  print 'Saving hdf5 file for fast gene loading...'
  save_hdf5_genes(E, gene_list, main_spring_dir + 'counts_norm_sparse_genes.hdf5')

  ##############
  print 'Saving hdf5 file for fast cell loading...'
  save_hdf5_cells(E,main_spring_dir + 'counts_norm_sparse_cells.hdf5')

  save_sparse_npz(E,main_spring_dir + 'counts_norm.npz', compressed = False)
  t0 = time.time()

  # subplot - original
  print 'making spring subplot - original ...'
  save_path = main_spring_dir + 'original'

  # make_spring_subplot() is the core function. feel free to change the params for your own purposes
  out = make_spring_subplot(E, gene_list, save_path,
                    normalize = True, tot_counts_final = total_counts,
                    min_counts = 1, min_cells = 3, min_vscore_pctl = 85,show_vscore_plot = False,
                    num_pc = 30,
                    k_neigh = 4,
                    num_force_iter = 100,
                    # the following is added for custom tsne and categorious by lizc@20180822
                    cell_groupings = cell_groupings,
                    custom_colors = {},
                    )
  
  print 'saving cell_filter.npy and cell_filter.txt ...'
  np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
  np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')
  print 'original subplot url: http://192.168.1.100:8000/springViewer_1_6_dev.html?datasets/' + os.path.basename(output_path.strip('/')) + '/original'
  
  # subplot custom
  if tsne_data:
    print 'making spring subplot - custom ...'
    save_path =main_spring_dir + 'custom'
  
    # make_spring_subplot() is the core function. feel free to change the params for your own purposes
    out = make_spring_subplot(E, gene_list, save_path,
                      normalize = True, tot_counts_final = total_counts,
                      min_counts = 1, min_cells = 3, min_vscore_pctl = 85,show_vscore_plot = False,
                      num_pc = 30,
                      k_neigh = 4,
                      num_force_iter = 100,
                      # the following is added for custom tsne and categorious by lizc@20180822
                      cell_groupings = cell_groupings,
                      custom_colors = {},
                      )
    # replace coordinates.txt with tsne_data
    
    print 'replacing coordinates.txt with tsne_data ...'
    positions = tsne_data['data']
    positions = positions * 5.0
    positions = positions - np.min(positions, axis = 0) - np.ptp(positions, axis = 0) / 2.0
    positions[:,0] = positions[:,0]  + 750
    positions[:,1] = positions[:,1]  + 250
    np.savetxt(save_path + '/coordinates.txt', np.hstack((np.arange(positions.shape[0])[:,None], positions)), fmt='%i,%.5f,%.5f')
    print 'saving cell_filter.npy and cell_filter.txt ...'
    np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
    np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')
    print 'custom subplot url: http://192.168.1.100:8000/springViewer_1_6_dev.html?datasets/' + os.path.basename(output_path.strip('/')) + '/custom'

  print 'Finished in %i seconds' %(time.time() - t0)
  

if __name__ == '__main__':
  parser = ArgumentParser(description = 'Generate dataset from 10X outs derived from "cellranger count", eg. outs/filtered_gene_bc_matrices/__YOUR_REFERENCE__/')
  parser.add_argument('input_path', help = '[required] path to 10X outs dir, eg /path/to/10x/outs/filtered_gene_bc_matrices/__YOUR_REFERENCE__/')
  parser.add_argument('output_path', help = '[required] path to SPRING datasets output directory, eg. /data2/lzc/SPRING_dev-master/datasets/__YOUR_PROJECT_NAME__')
  parser.add_argument('-m', action = 'store', dest = 'meta_data', help = '[optional] .csv file including cell groupings info, eg. annotation data from pbmc@meta.data. cellnames (that is rownames) intersected with cellranger results will be retained, colnames will be treated as group names')
  parser.add_argument('-t', action = 'store', dest = 'tsne_data', help = '[optional] .csv file including cell tSNE coordinate info, which will replace spring map, eg. dimension data from pbmc@dr$tsne@cell_loadings? cellnames (that is rownames) intersected with cellranger results will be retained, colnames will be treated as group names')
  
  args = parser.parse_args()
  if not (args.input_path and args.output_path):
    parser.print_help()
    exit(1)
  generate_data(args.input_path, args.output_path, args.meta_data, args.tsne_data)
