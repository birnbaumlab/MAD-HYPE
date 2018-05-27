import json

class SequencingData(object):
  def __init__(self, well_data=None, metadata=None, path=None):
    if path is not None:
      self.load_data(path)
    else:
      self.well_data = well_data
      self.metadata = metadata

  def load_data(self, path):
    data = json.load(open(path, 'r'))
    self.well_data = data['well_data']
    self.metadata = data['metadata']
    self.metadata['cells'] = [(tuple(alist), tuple(blist)) for alist, blist in self.metadata['cells']]
  def save_data(self, path):
    data = {'well_data': self.well_data, 'metadata': self.metadata}
    json.dump(data, open(path, 'w'))
  def save_data_R(self, alpha_path, beta_path):
    # saves data in a format usable with the R alphabetr implementation
    max_alpha = max([a for alist,_ in self.metadata['cells'] for a in alist])
    max_beta = max([b for _,blist in self.metadata['cells'] for b in blist])

    alpha_data,beta_data = [],[]
    for alist, blist in self.well_data:
      alpha_row, beta_row = [0]*(max_alpha+1), [0]*(max_beta+1)
      for a in alist:  alpha_row[a] = 1
      for b in blist:  beta_row[b] = 1
      alpha_data.append(alpha_row)
      beta_data.append(beta_row)

    f = open(alpha_path, 'w')
    f.write('\n'.join([','.join([str(v) for v in arow]) for arow in alpha_data]))
    f.close()

    f = open(beta_path, 'w')
    f.write('\n'.join([','.join([str(v) for v in brow]) for brow in beta_data]))
    f.close()


  def get_well_data(self, well_id = None):
    if well_id == None:
      return self.well_data
    else:
      return self.well_data[well_id]

  def get_metadata(self):
    return self.metadata
