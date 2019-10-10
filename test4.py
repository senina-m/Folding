# TensorFlow и tf.keras
import tensorflow as tf
from tensorflow import keras

# Вспомогательные библиотеки
import numpy as np
import matplotlib.pyplot as plt

print(tf.__version__)

fashion_mnist = keras.datasets.fashion_mnist
# in my case i need to download FASTA sequence
# and convert it to the one-hot code
fasta_arr = []

(train_images, train_labels), (test_images, test_labels) = fashion_mnist.load_data()

class_names = ['T-shirt/top', 'Trouser', 'Pullover', 'Dress', 'Coat',
               'Sandal', 'Shirt', 'Sneaker', 'Bag', 'Ankle boot']

train_images = train_images / 255.0
test_images = test_images / 255.0
# I don't need to convert and devide it

# here I have to know the length that my convolution layer can return
# TODO: give length_Of_window the right value
length_Of_window = 30
alphabet_len = 20
# 
model = keras.Sequential([
    keras.layers.Conv1D(input_shape=(len(fasta_arr), n, 1), )
    keras.layers.Dense(128, activation='relu'),
    keras.layers.Dense(10, activation='softmax')
])

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])

model.fit(train_images, train_labels, epochs=10)

test_loss, test_acc = model.evaluate(test_images, test_labels)

print('\nТочность на проверочных данных:', test_acc)

predictions = model.predict(test_images)

img = test_images[0]

img = (np.expand_dims(img,0))

predictions_single = model.predict(img)

print(predictions_single)
print(np.argmax(predictions_single[0]))
