import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat  #读取mat文件中数据的
import scipy.optimize as opt  #优化参数的，即训练使代价函数最小
from sklearn.metrics import classification_report  #用于评估分类器的好坏



#可视化数据
def load_mat(path):
    data = loadmat('ex4data1.mat')  # return a dict
    X = data['X']
    y = data['y'].flatten()
    return X, y


def plot_100_images(X):
    """随机画100个数字"""
    index = np.random.choice(range(5000), 100)
    images = X[index]
    fig, ax_array = plt.subplots(10, 10, sharey=True, sharex=True, figsize=(8, 8))  #10,10就是指一幅图里有100个小图
    for r in range(10):
        for c in range(10):
            ax_array[r, c].matshow(images[r*10 + c].reshape(20,20), cmap='gray_r')  #注意r和c的顺序
    plt.xticks([])
    plt.yticks([])
    plt.show()

X,y = load_mat('ex4data1.mat')
plot_100_images(X)



#读取数据
    #将y用0和1表示
def expand_y(y):  #将y用0和1表示，例如y[0]=6转化为y[0]=[0,0,0,0,0,1,0,0,0,0]。
    result = []
    for i in y:
        y_array = np.zeros(10)
        y_array[i-1] = 1
        result.append(y_array)
    return np.array(result)

    #获取训练数据
raw_X, raw_y = load_mat('ex4data1.mat')  
X = np.insert(raw_X, 0, 1, axis=1)  #插入的数据可以不用np.生成列，直接1简写就行
y = expand_y(raw_y)
X.shape, y.shape  #期待输出((5000, 401), (5000, 10))

    #展开参数  使用高级优化方法来优化神经网络时，我们需要将多个参数矩阵展开，再恢复
def serialize(a, b):
    return np.r_[a.flatten(),b.flatten()]  #np.r指连接两个展开的array

def deserialize(seq):  #再恢复参数。取决于几层神经网络，每层各有几个激活值
    return seq[:25*401].reshape(25, 401), seq[25*401:].reshape(10, 26)



#前馈（得到每层的输入输出）
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

    #获得每层的输入输出
def feed_forward(theta, X,):  #获得输入输出
    t1, t2 = deserialize(theta)
    a1 = X.T
    z2 = np.dot(t1,a1)
    a2 = np.insert(sigmoid(z2), 0, 1, axis=0)
    z3 = np.dot(t2,a2)
    a3 = sigmoid(z3)
    return a1, z2, a2, z3, a3

    #定义代价函数
def cost(theta, X, y):
    a1, z2, a2, z3, h = feed_forward(theta, X)  
    J=-y*np.log(h.T)-(1-y)*np.log(1-h.T)  #代价函数参考课程
    return J.sum()/len(X)

   #定义正则化后的代价函数
def regularized_cost(theta, X, y, l=1):
    t1, t2 = deserialize(theta)
    reg = np.sum(t1[:,1:] ** 2) + np.sum(t2[:,1:] ** 2)   #正则化的公式参考课程
    return l / (2 * len(X)) * reg + cost(theta, X, y)



#反向传播
    #定义S函数的导数
def sigmoid_gradient(z):
    return sigmoid(z) * (1 - sigmoid(z))

    #随机初始化
def random_init(size):  #总之维度不对就转置即可，自己对应好维度
    '''从范围中随机返回size个的值'''
    return np.random.uniform(-0.12, 0.12, size)

    #反向传播计算w的导数
def gradient(theta, X, y):
    t1, t2 = deserialize(theta)  #由于下面计算d2的时候要用到t2，所以要先定义一下
    a1, z2, a2, z3, h = feed_forward(theta, X)
    d3 = h.T - y  #(5000, 10)  #dn是指第n层的误差
    d2 = np.dot(d3,t2[:,1:] )* sigmoid_gradient(z2.T)  #(5000,25)
    D3=np.dot(a2,d3).T  # (10, 26)  #D3是指代价函数关于传到第3层的权重的偏导数/改变率。要不要转置，变成权重的维度一样就行
    D2=np.dot(a1,d2).T  # (25, 401)
    D = (1 / len(X)) * serialize(D2, D3)  # (10285,)
    return D



#正则化神经网络
def regularized_gradient(theta, X, y, l=1):
    """不惩罚偏置单元的参数"""
    D2, D3 = deserialize(gradient(theta, X, y)) #获得代价函数关于各个权重的导数
    t1,t2 = deserialize(theta)
    t1[:,0] = 0  #即把偏置单元的参数变成0
    t2[:,0] = 0
    reg_D2 = D2 + (l / len(X)) * t1
    reg_D3 = D3 + (l / len(X)) * t2
    return serialize(reg_D2, reg_D3)



#优化参数（即用现有的优化算法进行梯度下降）
def nn_training(X, y):
    init_theta = random_init(10285)  # 25*401 + 10*26
    res = opt.minimize(fun=regularized_cost,
                       x0=init_theta,
                       args=(X, y, 1),
                       method='TNC',
                       jac=regularized_gradient,
                       options={'maxiter': 400})
    return res



res = nn_training(X, y)  #大概优化20秒
res  #res.x是优化后的权重
    #检验优化后的效果/精确度
def accuracy(theta, X, y):
    _, _, _, _, h = feed_forward(res.x, X)
    y_pred = np.argmax(h.T, axis=1) + 1
    print(classification_report(y, y_pred))

accuracy(res.x, X, raw_y)  #precision和recall，一个是预测对的的确是对的概率，一个是对的被预测成对的概率，f1score是综合平均
