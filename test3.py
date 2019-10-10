import tensorflow as tf

# Создадим граф
graph = tf.get_default_graph()

# Создадим входное значение
input_value = tf.constant(1.0)

# Создадим переменную
weight = tf.Variable(0.8)

# Создаим выходное значение
output_value = weight * input_value

# создаим сессию
sess = tf.Session()

# Создаем оптимизиатор
desired_output_value = tf.constant(2.0)
loss = (output_value - desired_output_value) ** 2  # Функция ошибки
optim = tf.train.GradientDescentOptimizer(learning_rate=0.025)  # Оптимизатор
grads_and_vars = optim.compute_gradients(loss)

# Инициализируем переменные
init = tf.global_variables_initializer()
sess.run(init)

# Обучаем
train_step = tf.train.GradientDescentOptimizer(0.025).minimize(loss)
for i in range(100):
    sess.run(train_step)

print(sess.run(output_value))
