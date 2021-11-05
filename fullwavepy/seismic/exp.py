"""
Handle for various (meta)data associated 
with specific seismic experiments.

To add your own experiment, 
subclass the abstract `Experiment` class:
  class MyOwn(Experiment):
    pass
"""
from abc import ABC, abstractmethod
# from autologging import logged, traced


class Experiment(ABC):
  @abstractmethod
<<<<<<< HEAD
=======
  def _init_base_files(self):
    pass
  @abstractmethod
>>>>>>> 8b70a67d53c9b756ccbc6f04530d314d35991b08
  def _init_paths(self):
    pass
  @abstractmethod
  def _init_segy_mapping(self):
    pass
  @classmethod
  def create(cls, ID, *args, **kwargs):
    if ID == 'proteus':
      return PROTEUS(*args, **kwargs)
    elif ID is None:
      return None
    else:
      raise ValueError

