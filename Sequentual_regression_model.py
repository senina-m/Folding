import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


all_dataset = get_dataset()

def build_model():
    model = keras.Sequential([
        layers.Dense(20 * 50, activation='relu', input_shape=[len(train_dataset.keys())]),
        layers.Dense(20 * 50, activation='relu'),
        layers.Dense(50 * 50, activation='relu'),
        layers.Dense(50 * 50, activation='relu'),
        layers.Dense(50 * 50)
    ])

    optimizer = tf.keras.optimizers.RMSprop(0.001)

    model.compile(loss='mse',
                  optimizer=optimizer,
                  metrics=['mae', 'mse'])
    return model


model = build_model()
model.summary()


