import numpy as np 
import pandas as pd
import matplotlib as mplt
import matplotlib.pyplot as plt

class ModelUtils:

  @staticmethod
  def plot_model(x_train, x_test, y_train, y_test, model):
    """
    input model and actual data, plot error vs. y actual
    """
    fig, axes = plt.subplots(2)
    fig.dpi=100
    axes[0].scatter(y_train, model.predict(x_train).reshape(1,-1)-y_train, color = 'r',label='training')
    axes[1].scatter(y_test,model.predict(x_test).reshape(1,-1)-y_test, color = 'blue',label = 'test')
    fig.legend()
    plt.suptitle("error vs. y actual" )
    
  @staticmethod
  def combine_columns(nd_arrays):
    """
    combine nd array as columns
    nd_arrays: tuple of nd arrays of same dimension
    return: np.array
    """
    return np.column_stack((x_str,y_wl,Y_category))
    
  @staticmethod
  def subsampling(dataset, class_column_index, class_max_count, class_dict):
    """
    subsample each class to make sure every class has closer counts
    
    dataset -- nd.array: includes x column, y(category) etc
    class_column_index -- int: index of the category data in dataset
    class_max_count -- int: maximum number of counts allowed in each class 
    class_dict -- dict: {class: count}
    
    return subsampled data -- nd array
    """
    out = []
    for row in dataset:
        cls = row[class_column_index]
        rInt = np.random.randint(0, class_dict[cls])
        if rInt <= class_max_count:
            out.append(row)
    return np.array(out)
