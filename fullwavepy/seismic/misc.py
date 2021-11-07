from abc import ABC, abstractmethod

class BoxFactory:
  @classmethod
  def create(self, *args):
    if len(args) == 2:
      d = Box1d(*args)
    elif len(args) == 4:
      d = Box2d(*args)      
    elif len(args) == 6:
      d = Box3d(*args)   
    else:
      raise TypeError('Provide 2, 4 or 6 coordinates!')
    return d
class Box(ABC):
  """
  Spatial domain of a cuboidal shape.
  """  
  pass
class Box1d(Box):
  def __init__(self, x1, x2):
    self.box = [x1, x2]
    self.extent = [[x1, x2]]
    self.x1, self.x2 = self.box 
class Box2d(Box):
  def __init__(self, x1, x2, y1, y2):
    self.box = [x1, x2, y1, y2]
    self.extent = [[x1, x2], [y1, y2]]
    self.x1, self.x2, self.y1, self.y2 = self.box 
class Box3d(Box):
  def __init__(self, x1, x2, y1, y2, z1, z2):
    self.box = [x1, x2, y1, y2, z1, z2]
    self.extent = [[x1, x2], [y1, y2], [z1, z2]]
    self.xy = [[x1, x2], [y1, y2]]
    self.x1, self.x2, self.y1, self.y2, self.z1, self.z2 = self.box  
  def dims(self, dx):
    x1, x2, y1, y2, z1, z2 = self.box
    assert x2 >= x1
    assert y2 >= y1
    assert z2 >= z1
    nx1 = int((x2 - x1) / dx) + 1 
    nx2 = int((y2 - y1) / dx) + 1  
    nx3 = int((z2 - z1) / dx) + 1     
    return nx1, nx2, nx3
