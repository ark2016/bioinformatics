import pandas as pd
import numpy as np
from scipy.interpolate import BSpline, make_interp_spline
import string
from tqdm.notebook import tqdm

ALPHABET = string.ascii_uppercase
D = len(ALPHABET) # Размерность пространства
SPLINE_DEGREE = 3
INTEGRATION_POINTS = 100 
D_MAX_ESTIMATE = 100 # подобрать
BETA_LENGTH_FACTOR = 5.0 # подобрать

char_to_onehot = {char: np.zeros(D) for char in ALPHABET}
for i, char in enumerate(ALPHABET):
    char_to_onehot[char][i] = 1

def encode_sequence(sequence):
    """Преобразует строку в последовательность 26-мерных one-hot векторов."""
    encoded = []
    for char in sequence:
        if char in char_to_onehot:
            encoded.append(char_to_onehot[char])
        else:
            encoded.append(np.zeros(D))
    return np.array(encoded)

def create_bspline_for_sequence(encoded_sequence, degree=SPLINE_DEGREE):
    """
    Создает 26 одномерных B-сплайнов для закодированной последовательности.
    Возвращает список из 26 объектов BSpline.
    """
    L = encoded_sequence.shape[0]
    if L < degree + 1:
        if L == 0: return None
        if L == 1:
             return [BSpline.construct_fast([0,1], [encoded_sequence[0,j], encoded_sequence[0,j]], 0) for j in range(D)]

    t_param = np.linspace(0, 1, L)

    spline_components = []
    for j in range(D):
        y_coords = encoded_sequence[:, j]
        spl = make_interp_spline(t_param, y_coords, k=degree)
        spline_components.append(spl)
    return spline_components

def compute_l2_distance_between_spline_parts(spline_short_components, spline_long_sub_components, integration_points=INTEGRATION_POINTS):
    """
    Вычисляет интегральное L2 расстояние между двумя наборами сплайнов.
    Использует численное интегрирование.
    """
    if spline_short_components is None or spline_long_sub_components is None:
        return np.inf
    
    t_eval = np.linspace(0, 1, integration_points)
    total_squared_diff = 0.0

    for j in range(D):
        val_short = spline_short_components[j](t_eval)
        val_long_sub = spline_long_sub_components[j](t_eval)
        total_squared_diff += np.sum((val_short - val_long_sub)**2)

    l2_distance = np.sqrt(total_squared_diff / integration_points)
    return l2_distance

def create_scaled_sub_spline_components(full_long_spline_components, ts_start, lambda_ratio):
    """
    Создает масштабированные компоненты сплайна для подинтервала.
    Возвращает список из 26 объектов BSpline, представляющих часть
    длинного сплайна, масштабированную на [0, 1].
    """
    ts_end = ts_start + lambda_ratio
    scaled_components = []
    for j in range(D):
        original_spl = full_long_spline_components[j]
        num_eval_points = 50
        t_sub_original = np.linspace(ts_start, ts_end, num_eval_points)
        y_sub_original = original_spl(t_sub_original)
        t_scaled = np.linspace(0, 1, num_eval_points)
        degree = min(SPLINE_DEGREE, num_eval_points - 1)
        if degree < 1:
            scaled_components.append(BSpline.construct_fast([0,1], [y_sub_original[0], y_sub_original[0]], 0))
        else:
            scaled_components.append(make_interp_spline(t_scaled, y_sub_original, k=degree))
    return scaled_components

def compare_strings_with_splines(df_sequences, template_index, long_string_index):
    """
    Сравнивает сплайн-шаблон с частями длинного сплайна.
    Возвращает минимальное L2 расстояние, M_form, M_length, M_total.
    """
    seq_template = df_sequences.iloc[template_index]
    seq_long = df_sequences.iloc[long_string_index]

    L_short = len(seq_template)
    L_long = len(seq_long)

    if L_short == 0 or L_long == 0 or L_short > L_long:
        return np.inf, 0.0, 0.0, 0.0
    
    encoded_template = encode_sequence(seq_template)
    encoded_long = encode_sequence(seq_long)

    spline_template_components = create_bspline_for_sequence(encoded_template)
    spline_long_components = create_bspline_for_sequence(encoded_long)

    if spline_template_components is None or spline_long_components is None:
        return np.inf, 0.0, 0.0, 0.0

    lambda_ratio = L_short / L_long
    min_distance = np.inf
    num_shifts = L_long - L_short + 1

    if num_shifts <= 0:
        ts_start_points = [0.0]
    else:
        ts_start_points = np.linspace(0, 1 - lambda_ratio, num_shifts)

    for ts_start in ts_start_points:
        spline_long_sub_components = create_scaled_sub_spline_components(
            spline_long_components, ts_start, lambda_ratio
        )
        current_distance = compute_l2_distance_between_spline_parts(
            spline_template_components, spline_long_sub_components
        )
        if current_distance < min_distance:
            min_distance = current_distance

    M_form = 1 - (min_distance / D_MAX_ESTIMATE)
    M_form = max(0.0, min(1.0, M_form))
    
    M_length = np.exp(-BETA_LENGTH_FACTOR * np.abs(lambda_ratio - 1))

    M_total = M_form * M_length

    return min_distance, M_form, M_length, M_total
