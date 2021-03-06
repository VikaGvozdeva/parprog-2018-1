﻿# Каталог студента Гвоздевой Виктории Евгеньевны

## Структура каталога

- 1-test-version - каталог для размещения файлов тестовой версии
- 2-openmp - каталог для размещения файлов OpenMP-версии
- 3-tbb - каталог для размещения файлов TBB-версии

Умножение разреженных матриц. Элементы типа double. Формат хранения матрицы – строковый (CRS).
Задача: A * B = C
Алгоритм умножения матриц известен. Строки матрицы A умножаются на столбцы матрицы B. В результате получаем матрицу C.

Для хранения матрицы в формате CRS используются три массива. Первый (val) массив хранит значения элементов построчно 
(строки рассматриваются по порядку сверху вниз), второй (col_ind) – номера столбцов для каждого элемента, а третий (row_ptr) заменяет номера
строк, используемые в координатном формате, на индекс начала каждой строки. Количество элементов массива row_ptr равно N + 1, где N - размерность
матрицы. При этом элементы строки i в массиве val находятся по индексам от row_ptr[i] до row_ptr[i + 1] – 1 включительно (i-ый элемент
массива row_ptr указывает на начало i-ой строки).
Данная реализация выполняет следующие действия для умножения матриц:
1) Матрица B транспонируется, так как в данном формате доступ к строке осуществляется проще, чем к столбцу.
2) Так как при умножение двух разреженных матриц результирующая матрица может получиться как разреженной, так и плотной, сначала получаем 
портрет итоговой матрицы.
3) Умножение двух матриц с уже известным портретом результирующей матрицы.
