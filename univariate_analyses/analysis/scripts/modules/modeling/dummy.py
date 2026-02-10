import numpy as np

class BaselineModel:
    '''
    Fit and predict methods for a unconditional mean model
    '''

    def __init__(self, **kwargs):

        self.mean_ = None

    def fit(self, X, y):

        self.mean_ = np.mean(y)

    def predict(self, X):
        
        return np.full(shape=(X.shape[0],), fill_value=self.mean_)
