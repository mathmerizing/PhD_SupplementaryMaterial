import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import imageio

# Source 1: https://github.com/dynamicslab/databook_python/blob/master/CH01/CH01_SEC06_1.ipynb
# Source 2: https://github.com/dynamicslab/databook_python/blob/master/CH01/CH01_SEC06_2_3_4.ipynb
# NOTE: This code is heavily based on the supplementary code to the book "Data-driven science and engineering: Machine learning, dynamical systems, and control" by Brunton and Kutz (2022)

# if file allFaces.mat does not exist:
#   download the dataset from https://github.com/dynamicslab/databook_python/raw/master/DATA/allFaces.mat
#   and save it in the same folder as this file

# download the dataset
if not os.path.isfile('allFaces.mat'):
    import urllib.request
    url = 'https://github.com/dynamicslab/databook_python/raw/master/DATA/allFaces.mat'
    urllib.request.urlretrieve(url, 'allFaces.mat')
    print('Downloaded the dataset')

plt.rcParams['figure.figsize'] = [10, 10]
plt.rcParams.update({'font.size': 18})

mat_contents = scipy.io.loadmat('allFaces.mat')
faces = mat_contents['faces']
m = int(mat_contents['m'])
n = int(mat_contents['n'])
nfaces = np.ndarray.flatten(mat_contents['nfaces'])

# allPersons = np.zeros((n*6,m*6))
# count = 0

# for j in range(6):
#     for k in range(6):
#         allPersons[j*n : (j+1)*n, k*m : (k+1)*m] = np.reshape(faces[:,np.sum(nfaces[:count])],(m,n)).T
#         count += 1
        
# img = plt.imshow(allPersons)
# img.set_cmap('gray')
# plt.axis('off')
# plt.show()

# We use the first 36 people for training data
trainingFaces = faces[:,:np.sum(nfaces[:36])]
avgFace = np.mean(trainingFaces,axis=1) # size n*m by 1

# Compute eigenfaces on mean-subtracted training data
X = trainingFaces - np.tile(avgFace,(trainingFaces.shape[1],1)).T
U, S, VT = np.linalg.svd(X,full_matrices=0)

# plot the eigenvalues decay
plt.semilogy(S,'-o',linewidth=2)
plt.xlabel('Component number')
plt.ylabel('Singular value')
plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(121)
img_avg = ax1.imshow(np.reshape(avgFace,(m,n)).T)
img_avg.set_cmap('gray')
plt.axis('off')

ax2 = fig1.add_subplot(122)
img_u1 = ax2.imshow(np.reshape(U[:,0],(m,n)).T)
img_u1.set_cmap('gray')
plt.axis('off')

plt.show()



## Now show eigenface reconstruction of image that was omitted from test set
testFace = faces[:,np.sum(nfaces[:36])] # First face of person 37
plt.imshow(np.reshape(testFace,(m,n)).T)
plt.set_cmap('gray')
plt.title('Original Image')
plt.axis('off')
plt.show()

# load julian_face.jpg into a numpy array
julian_face = imageio.imread('julian_face.jpg')
julian_face = julian_face[:,:,0].T.flatten()
print(julian_face.shape)
print(testFace.shape)
img = plt.imshow(np.reshape(julian_face,(m,n)).T)
img.set_cmap('gray')
plt.axis('off')
plt.show()
plt.show()

# # load thomas_face.jpg into a numpy array
# thomas_face = imageio.imread('thomas_face.jpg')
# thomas_face = thomas_face.T.flatten()
# print(thomas_face.shape)
# print(testFace.shape)
# img = plt.imshow(np.reshape(thomas_face,(m,n)).T)
# img.set_cmap('gray')
# plt.axis('off')
# plt.show()
# plt.show()

# testFaceMS = testFace - avgFace
# r_list = [25, 50, 100, 200, 400, 800, 1600]

# for r in r_list:
#     reconFace = avgFace + U[:,:r]  @ U[:,:r].T @ testFaceMS
#     img = plt.imshow(np.reshape(reconFace,(m,n)).T)
#     img.set_cmap('gray')
#     plt.title('r = ' + str(r))
#     plt.axis('off')
#     plt.show()

test_face_MS = julian_face - avgFace
r_list = [25, 50, 100, 200, 400, 800, 1600]

for r in r_list:
    reconFace = avgFace + U[:,:r]  @ U[:,:r].T @ test_face_MS
    img = plt.imshow(np.reshape(reconFace,(m,n)).T)
    img.set_cmap('gray')
    plt.title('r = ' + str(r))
    plt.axis('off')
    plt.show()

