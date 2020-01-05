import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# see built_data.py
all_dataset = get_dataset()
train_data, test_data = split_data(all_dataset)


def build_model():
    model = keras.Sequential([
        layers.Dense(20 * 50, activation='relu', input_shape=[len(train_dataset.keys())]),
        layers.Dense(100 * 50, activation='relu'),
        layers.Dense(70 * 50, activation='relu'),
        layers.Dense(50 * 50, activation='relu'),
        layers.Dense(50 * 50, activation='linear')
    ])

    optimizer = tf.keras.optimizers.RMSprop(0.001)

    model.compile(loss='mse',
                  optimizer=optimizer,
                  metrics=['mae', 'mse'])
    return model


model = build_model()
model.summary()
example_batch = train_data[:10]
example_result = model.predict(example_batch),


