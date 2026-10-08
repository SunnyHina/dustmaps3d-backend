"""Scientific algorithms ported from the original dustmap/function.py.

Input datasets and output locations are explicit; no web-framework dependency.
"""
import os
import time
import uuid
import shutil
import tempfile
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy import units as u
from astropy_healpix import HEALPix
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import LogLocator, ScalarFormatter, FuncFormatter, MultipleLocator, FixedLocator
from matplotlib.patches import Circle, Rectangle
import matplotlib.patheffects as pe
import cmasher as cmr
import warnings
import matplotlib.colors as mcolors
import io
import base64
from astropy.coordinates import SkyCoord, CartesianRepresentation
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore", message="overflow encountered in exp")

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("4_colormap", ["green",  "blue", "yellow", "red"])
cmap4 = LinearSegmentedColormap.from_list("4_colormap", ["blue",  "green", "yellow", "red"])
cmap5 = LinearSegmentedColormap.from_list("4_colormap", ["blue", "green", "tan" , "yellow", "red"])


DPR_EPS = 1e-12
DPR_CYLINDER_HALF_HEIGHT_KPC = 0.05
DPR_DEFAULT_EXTRA_CAP_SECTION_IDS = {7}
DPR_SB21_CAP_SECTION_IDS = {21}
DPR_SB_LABEL_FONTSIZE = 6.4


def bubble_diffuse(x, h, b_lim, diffuse_dust_rho, bubble):
    """
    气泡和弥散介质函数
    确保所有参数都是标量值
    """
    
    span = 0.01
    span_0 = h / np.sin(np.deg2rad(np.abs(b_lim)))
    Cum_EBV_0 = span_0 * diffuse_dust_rho
    C_0 = Cum_EBV_0 * (1 - np.exp(- (bubble) / span_0))
    f = (Cum_EBV_0 * (1 - np.exp(-x / span_0))) - C_0
    a = 1/np.exp(5 * bubble /span)
    b = 1 / (1 + np.exp(5 * bubble/span))
    c = 0.5
    deta = C_0/((1+a)*(c-b))
    return np.where(x < (bubble), 0, f) + deta*(1+a)*((1 / (1 + np.exp(-5 * ((x - bubble)/span))) )-b)


def component4(x, b_lim, bubble, diffuse_dust_rho, h, distance_1, span_1, Cum_EBV_1, 
               distance_2, span_2, Cum_EBV_2, distance_3, span_3, Cum_EBV_3, 
               distance_4, span_4, Cum_EBV_4):
    """
    4组件消光模型
    """
    
    Numerator_1 = Cum_EBV_1*(1/np.exp(5 * (distance_1 + (span_1*2) + bubble) /span_1) + 1)
    Numerator_2 = Cum_EBV_2*(1/np.exp(5 * (distance_2 + (span_2*2) + bubble)/span_2) + 1)
    Numerator_3 = Cum_EBV_3*(1/np.exp(5 * (distance_3 + (span_3*2) + bubble)/span_3) + 1)
    Numerator_4 = Cum_EBV_4*(1/np.exp(5 * (distance_4 + (span_4*2) + bubble)/span_4) + 1)
    
    return (bubble_diffuse(x, h, b_lim, diffuse_dust_rho, bubble)
                     
            +((Numerator_1/ (1 + np.exp(-5 * ((x) - (distance_1 + (span_1*2) + bubble))/span_1))) 
            -(Numerator_1 / (1 + np.exp(5 * (distance_1 + (span_1*2) + bubble)/span_1))))
            
            +((Numerator_2 / (1 + np.exp(-5 * ((x) - (distance_2 + (span_2*2) + bubble))/span_2))) 
            -(Numerator_2 / (1 + np.exp(5 * ((distance_2 + (span_2*2) + bubble))/span_2))))

            +((Numerator_3 / (1 + np.exp(-5 * ((x) - (distance_3 + (span_3*2) + bubble))/span_3))) 
            -(Numerator_3 / (1 + np.exp(5 * ((distance_3 + (span_3*2) + bubble))/span_3))))

            +((Numerator_4 / (1 + np.exp(-5 * ((x) - (distance_4 + (span_4*2) + bubble))/span_4))) 
            -(Numerator_4 / (1 + np.exp(5 * ((distance_4 + (span_4*2) + bubble))/span_4))))
            )
 
def diffusion_derived_function(x, b_lim, diffuse_dust_rho, h):
    """
    弥散函数的导数
    """
    
    span_0 = h / np.sin(np.deg2rad(np.abs(b_lim)))
    return diffuse_dust_rho * (np.exp(- x / span_0))

def sigmoid(x, a, b, c):
    return c / (1 + np.exp(-b * (x - a)))

def derivative_of_sigmoid(x, a, b, c):
    if c == 0:
        return np.zeros_like(x)
    else:
        return b * c * sigmoid(x, a, b, 1) * (1 - (sigmoid(x, a, b, 1)))

def sigmoid_of_component(bubble, distance, span, Cum_EBV):
    a = distance + (2*span) + bubble
    b = 5 / span
    c = Cum_EBV*(1/np.exp(5 * a /span) + 1)
    return a, b, c

def derivative_of_component4(x, b_lim, bubble, diffuse_dust_rho, h, distance_1, span_1, Cum_EBV_1, 
                           distance_2, span_2, Cum_EBV_2, distance_3, span_3, Cum_EBV_3, 
                           distance_4, span_4, Cum_EBV_4):
    """
    4组件模型的导数
    """
    
    a_1, b_1, c_1 = sigmoid_of_component(bubble, distance_1, span_1, Cum_EBV_1)
    a_2, b_2, c_2 = sigmoid_of_component(bubble, distance_2, span_2, Cum_EBV_2)
    a_3, b_3, c_3 = sigmoid_of_component(bubble, distance_3, span_3, Cum_EBV_3)
    a_4, b_4, c_4 = sigmoid_of_component(bubble, distance_4, span_4, Cum_EBV_4)
    
    return (np.where(x < bubble, 0, diffusion_derived_function(x, b_lim, diffuse_dust_rho, h)) 
            + derivative_of_sigmoid(x, a_1, b_1, c_1) 
            + derivative_of_sigmoid(x, a_2, b_2, c_2) 
            + derivative_of_sigmoid(x, a_3, b_3, c_3) 
            + derivative_of_sigmoid(x, a_4, b_4, c_4))

ALLOWED_EXTENSIONS = {'csv', 'fits', 'fit'}

def allowed_file(filename):
    """检查文件扩展名"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def validate_and_prepare_dataframe(df):
    """
    自动检测坐标系，验证数据，并进行必要的转换。
    返回: (DataFrame, str_error_message)
    如果成功，返回处理后的DataFrame和None。
    如果失败，返回None和错误信息。
    """
    # 清理列名：转为小写并去除首尾空格，以便匹配
    df.columns = [col.lower().strip() for col in df.columns]
    if df.columns.duplicated().any():
        return None, "Duplicate column names are not supported."
    for column in ('l', 'b', 'ra', 'dec', 'd'):
        if column in df:
            try:
                df[column] = pd.to_numeric(df[column], errors='raise')
            except (ValueError, TypeError):
                return None, f"Column {column} must contain numeric values."
            if np.isinf(df[column]).any():
                return None, f"Column {column} must contain finite values."
    df_columns = set(df.columns)
    
    # 检查必需的距离列
    if 'd' not in df_columns:
        return None, "文件中必须包含距离列 'd'。"


    # 1. 自动检测坐标系
    if 'l' in df_columns and 'b' in df_columns:
        
        # 验证范围
        if not df['l'].between(0, 360).all():
            return None, "银经 (l) 的值必须在 0-360 度范围内。"
        if not df['b'].between(-90, 90).all():
            return None, "银纬 (b) 的值必须在 -90 到 90 度范围内。"

    elif 'ra' in df_columns and 'dec' in df_columns:
        
        # 验证范围
        if not df['ra'].between(0, 360).all():
            return None, "赤经 (ra) 的值必须在 0-360 度范围内。"
        if not df['dec'].between(-90, 90).all():
            return None, "赤纬 (dec) 的值必须在 -90 到 90 度范围内。"
        
        # 进行坐标转换
        try:
            coords = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg, frame='icrs')
            galactic_coords = coords.galactic
            df['l'] = galactic_coords.l.degree
            df['b'] = galactic_coords.b.degree
        except Exception as e:
            return None, f"赤道坐标转换失败: {str(e)}"

    else:
        return None, "文件列不符合要求。必须包含 ('l', 'b', 'd') 或 ('ra', 'dec', 'd')。"

    # 验证距离
    if not ((df['d'] > 0) | (df['d'].isna())).all():
        return None, "距离 (d) 的值必须是大于 0 的数字，或将该单元格留空。不支持 0 或负值。"
        
    # 确保返回的 DataFrame 包含 l, b, d
    return df, None

def calculate_dust_properties(l_coords, b_coords, distances, df_global):
    """
    统一的计算函数，适配处理后的数据 (data.parquet)。
    它接收银道坐标和距离，返回一个包含所有物理量的字典。
    如果输入的 distance 为空(NaN)，则使用该方向上的最大可靠距离进行计算。

    Args:
        l_coords (array-like): 银经 (Galactic Longitude l) in degrees.
        b_coords (array-like): 银纬 (Galactic Latitude b) in degrees.
        distances (array-like): 距离 (distance d) in kpc. 可能包含 NaN 值。
        df_global (pd.DataFrame): 预加载的处理后数据，索引为 pix_1024.

    Returns:
        dict: 包含计算结果列表的字典。
    """
    if df_global is None or df_global.empty:
        raise ValueError("Global dust map data is not loaded.")

    # 确保输入是 numpy 数组
    l_coords = np.atleast_1d(l_coords)
    b_coords = np.atleast_1d(b_coords)
    distances = np.atleast_1d(distances)
    
    nside = 1024
    # 将 l, b 坐标转换为 HEALPix 像素 ID
    pix_ids = hp.ang2pix(nside, l_coords, b_coords, lonlat=True, nest=False)

    # 准备结果容器
    results = {
        'E(B-V)': [],
        'dust_density': [],
        'sigma': [],
        'max_distance': []
    }

    # 直接使用 pix_id (即索引) 进行查找
    # .reindex() 可以优雅地处理存在和不存在的 pix_id
    # 找到的行将保留，未找到的行将填充为 NaN
    matched_rows = df_global.reindex(pix_ids)

    # 遍历每一行进行计算
    for (_, data_row), d_value in zip(matched_rows.iterrows(), distances):
        # 如果 data_row 是 NaN，说明 pix_id 不在数据中（虽然现在不太可能，因为数据是密集的）
        if data_row.isna().all():
            results['E(B-V)'].append(np.nan)
            results['dust_density'].append(np.nan)
            results['sigma'].append(np.nan)
            results['max_distance'].append(np.nan)
            continue
        
        max_dist = data_row['max_distance']
        
        d_for_calc = max_dist if pd.isna(d_value) else d_value

        ebv = component4(d_for_calc, data_row['b_lim'], data_row['bubble'], data_row['diffuse_dust_rho'], 
                         data_row['h'], data_row['distance_1'], data_row['span_1'], data_row['Cum_EBV_1'],
                         data_row['distance_2'], data_row['span_2'], data_row['Cum_EBV_2'],
                         data_row['distance_3'], data_row['span_3'], data_row['Cum_EBV_3'],
                         data_row['distance_4'], data_row['span_4'], data_row['Cum_EBV_4'])
        
        dust_density = derivative_of_component4(d_for_calc, data_row['b_lim'], data_row['bubble'], data_row['diffuse_dust_rho'], 
                                                data_row['h'], data_row['distance_1'], data_row['span_1'], data_row['Cum_EBV_1'],
                                                data_row['distance_2'], data_row['span_2'], data_row['Cum_EBV_2'],
                                                data_row['distance_3'], data_row['span_3'], data_row['Cum_EBV_3'],
                                                data_row['distance_4'], data_row['span_4'], data_row['Cum_EBV_4'])

        sigma_final = data_row['sigma']
        
        results['E(B-V)'].append(ebv)
        results['dust_density'].append(dust_density) # 单位 mag/kpc, 数值上等于 mmag/pc
        results['sigma'].append(sigma_final)
        results['max_distance'].append(max_dist)
        
    return results

def query_polygon_robust_helper(lon_deg, lat_deg, is_galactic_coord):
    """
    一个辅助函数，将经纬度顶点转换为query_polygon所需的三维向量。
    """
    if not is_galactic_coord:
        # 如果是赤道坐标，先转为银道坐标
        coords_icrs = SkyCoord(ra=lon_deg*u.deg, dec=lat_deg*u.deg, frame='icrs')
        coords_galactic = coords_icrs.transform_to('galactic')
        lon_rad = coords_galactic.l.radian
        lat_rad = coords_galactic.b.radian
    else:
        lon_rad = np.radians(lon_deg)
        lat_rad = np.radians(lat_deg)
    
    theta_rad = np.pi/2. - lat_rad
    phi_rad = lon_rad
    return hp.ang2vec(theta_rad, phi_rad)

def query_polygon_robust(nside, vertices_deg_lon, vertices_deg_lat, is_galactic):
    """
    一个稳健的多边形查询函数，能处理跨越0/360度边界的情况。
    """
    
    # 将顶点经纬度转换为三维向量的辅助函数
    def get_vec_from_lonlat(lon_deg, lat_deg, is_galactic_coord):
        if not is_galactic_coord:
            # 如果是赤道坐标，先转为银道坐标
            coords_icrs = SkyCoord(ra=lon_deg*u.deg, dec=lat_deg*u.deg, frame='icrs')
            coords_galactic = coords_icrs.transform_to('galactic')
            lon_rad = coords_galactic.l.radian
            lat_rad = coords_galactic.b.radian
        else:
            lon_rad = np.radians(lon_deg)
            lat_rad = np.radians(lat_deg)
        
        theta_rad = np.pi/2. - lat_rad
        phi_rad = lon_rad
        return hp.ang2vec(theta_rad, phi_rad)

    lon_min, lon_max = min(vertices_deg_lon), max(vertices_deg_lon)

    # 检查是否跨越了0/360度边界
    # (这里我们假设输入已经是处理过的，比如RA从-10变成350)
    # 比如用户输入 lon_min=350, lon_max=10
    # 我们需要检查原始输入，而不是这里的min/max
    # 一个简单的检查是，如果 lon_min_orig > lon_max_orig
    
    # *** 为了简化，我们假设函数调用者已经处理了跨界问题 ***
    # *** 这里的逻辑是，如果经度范围非常大，就可能跨界 ***
    # 一个更简单的处理方式是检查输入范围
    if vertices_deg_lon[0] > vertices_deg_lon[1]: # e.g., from 350 to 10
        # 跨界查询，拆分成两个多边形
        
        # 多边形1: [lon_min, 360] x [lat_min, lat_max]
        lon1 = np.array([vertices_deg_lon[0], 360, 360, vertices_deg_lon[0]])
        lat1 = vertices_deg_lat
        vec1 = get_vec_from_lonlat(lon1, lat1, is_galactic)
        pix1 = hp.query_polygon(nside, vec1, inclusive=True)
        
        # 多边形2: [0, lon_max] x [lat_min, lat_max]
        lon2 = np.array([0, vertices_deg_lon[1], vertices_deg_lon[1], 0])
        lat2 = vertices_deg_lat
        vec2 = get_vec_from_lonlat(lon2, lat2, is_galactic)
        pix2 = hp.query_polygon(nside, vec2, inclusive=True)
        
        # 合并两个区域的像素
        return np.union1d(pix1, pix2)
    else:
        # 正常查询
        vec = get_vec_from_lonlat(vertices_deg_lon, vertices_deg_lat, is_galactic)
        return hp.query_polygon(nside, vec, inclusive=True)

def hpix_Delta(df, distance_low, distance_up):
    """
    计算给定距离范围内的平均尘埃密度。
    OPTIMIZED VECTORIZED VERSION.

    Args:
        df (pd.DataFrame): 包含所需数据点的DataFrame，索引为 pix_1024。
        distance_low (float): 距离范围的下限 (kpc).
        distance_up (float): 距离范围的上限 (kpc).

    Returns:
        tuple: (pix_ids, dust_density_values)
    """
    if df.empty:
        return np.array([]), np.array([])

    # --- VECTORIZED CALCULATION ---
    # Instead of df.apply, we pass entire columns (Pandas Series) to the function.
    # `distance_low` and `distance_up` are scalars, NumPy handles broadcasting correctly.
    
    ebv_low = component4(
        distance_low, df['b_lim'], df['bubble'], df['diffuse_dust_rho'], 
        df['h'], df['distance_1'], df['span_1'], df['Cum_EBV_1'],
        df['distance_2'], df['span_2'], df['Cum_EBV_2'],
        df['distance_3'], df['span_3'], df['Cum_EBV_3'],
        df['distance_4'], df['span_4'], df['Cum_EBV_4']
    )

    ebv_up = component4(
        distance_up, df['b_lim'], df['bubble'], df['diffuse_dust_rho'], 
        df['h'], df['distance_1'], df['span_1'], df['Cum_EBV_1'],
        df['distance_2'], df['span_2'], df['Cum_EBV_2'],
        df['distance_3'], df['span_3'], df['Cum_EBV_3'],
        df['distance_4'], df['span_4'], df['Cum_EBV_4']
    )
    # --- END VECTORIZED CALCULATION ---

    delta_ebv = ebv_up - ebv_low
    
    delta_distance = distance_up - distance_low
    if delta_distance == 0:
        delta_distance = 1e-9

    dust_density = delta_ebv / delta_distance

    distance_midpoint = (distance_up + distance_low) / 2
    
    # `np.where` is also vectorized. It takes a boolean Series as a condition.
    final_dust_values = np.where(
        distance_midpoint < df['max_distance'],
        dust_density,
        hp.UNSEEN
    )
    
    return df.index.to_numpy(), final_dust_values

def get_pixels_in_galactic_box(nside, l_min, l_max, b_min, b_max):
    """
    【推荐】精确获取银道坐标系矩形区域内的所有HEALPix像素。
    采用全天像素反向查询法，逻辑简单，结果可靠，避免所有几何变形问题。

    Args:
        nside (int): HEALPix的nside参数。
        l_min (float): 银经下限 (degrees)。
        l_max (float): 银经上限 (degrees)。
        b_min (float): 银纬下限 (degrees)。
        b_max (float): 银纬上限 (degrees)。

    Returns:
        np.ndarray: 落在指定区域内的像素ID数组。
    """
    npix = hp.nside2npix(nside)
    all_pixels = np.arange(npix)

    theta, phi = hp.pix2ang(nside, all_pixels, nest=False)

    pix_l = np.degrees(phi)
    pix_b = 90.0 - np.degrees(theta)

    # 筛选银纬
    b_mask = (pix_b >= b_min) & (pix_b <= b_max)

    # 筛选银经，需要特别处理跨0/360度的情况
    if l_min > l_max:
        # 例如 l from 350 to 10. 逻辑是 (l >= 350) OR (l <= 10)
        l_mask = (pix_l >= l_min) | (pix_l <= l_max)
    else:
        # 例如 l from 10 to 350. 逻辑是 (l >= 10) AND (l <= 350)
        l_mask = (pix_l >= l_min) & (pix_l <= l_max)

    # 合并两个筛选条件
    final_mask = b_mask & l_mask
    
    # 返回通过所有筛选条件的像素ID
    return all_pixels[final_mask]

def get_pixels_in_equatorial_box(nside, ra_min, ra_max, dec_min, dec_max):
    """
    精确获取赤道坐标系矩形区域内的所有HEALPix像素。
    采用全天像素反向查询法，逻辑最简单，结果最可靠，避免所有几何变形问题。

    Args:
        nside (int): HEALPix的nside参数。
        ra_min (float): 赤经下限 (degrees)。
        ra_max (float): 赤经上限 (degrees)。
        dec_min (float): 赤纬下限 (degrees。
        dec_max (float): 赤纬上限 (degrees)。

    Returns:
        np.ndarray: 落在指定区域内的像素ID数组。
    """
    npix = hp.nside2npix(nside)
    all_pixels = np.arange(npix)

    # healpy返回的是 (theta, phi)，其中 theta 是极角 (colatitude)
    theta, phi = hp.pix2ang(nside, all_pixels, nest=False)
    
    pix_coords_gal = SkyCoord(l=(phi*u.rad), b=(np.pi/2. - theta)*u.rad, frame='galactic')
    pix_coords_icrs = pix_coords_gal.transform_to('icrs')
    pix_ra = pix_coords_icrs.ra.deg
    pix_dec = pix_coords_icrs.dec.deg

    # 筛选赤纬
    dec_mask = (pix_dec >= dec_min) & (pix_dec <= dec_max)

    # 筛选赤经，需要特别处理跨0/360度的情况
    if ra_min > ra_max:
        # 例如 RA from 350 to 10. 逻辑是 (RA >= 350) OR (RA <= 10)
        ra_mask = (pix_ra >= ra_min) | (pix_ra <= ra_max)
    else:
        # 例如 RA from 10 to 350. 逻辑是 (RA >= 10) AND (RA <= 350)
        ra_mask = (pix_ra >= ra_min) & (pix_ra <= ra_max)

    # 合并两个筛选条件
    final_mask = dec_mask & ra_mask
    
    # 返回通过所有筛选条件的像素ID
    return all_pixels[final_mask]

def plot_partial_dust_map_ic_improved(params, df_global, output_dir, public_base_url):
    """
    改进的绘图函数，解决赤道坐标系的各种问题
    """
    # 参数提取
    coord_system = params['coord_system']
    c1_min_orig, c1_max_orig = params['lon_min'], params['lon_max']
    c2_min_orig, c2_max_orig = params['lat_min'], params['lat_max']
    nside = 1024
    
    # 边界处理
    epsilon = 1e-6
    c2_max = min(c2_max_orig, 90.0 - epsilon)
    c2_min = max(c2_min_orig, -90.0 + epsilon)
    c1_max = min(c1_max_orig, 360.0 - epsilon)
    c1_min = max(c1_min_orig, 0.0 + epsilon)
    if c2_min >= c2_max:
        c2_min = c2_max - epsilon
    if c1_min >= c1_max:
        c1_min = c1_max - epsilon
    
    # 检查是否为全天查询
    lon_span = abs(c1_max_orig - c1_min_orig)
    lat_span = abs(c2_max_orig - c2_min_orig)
    is_full_sky = (lon_span >= 360.0) and (lat_span >= 179.0)
    
    # 获取像素ID
    if is_full_sky:
        # 全天查询
        pix_ids_in_range = np.arange(hp.nside2npix(nside))
    elif coord_system == 'galactic':
        c1_min = c1_min_orig % 360
        c1_max = c1_max_orig % 360

        if lon_span >= 360.0:
            theta_min_rad = np.pi/2. - np.radians(c2_max)
            theta_max_rad = np.pi/2. - np.radians(c2_min)
            pix_ids_in_range = hp.query_strip(nside, theta_min_rad, theta_max_rad, nest=False)
        else:
            # 使用更稳健的像素查询方法，避免query_polygon的退化角问题
            pix_ids_in_range = get_pixels_in_galactic_box(nside, c1_min, c1_max, c2_min, c2_max)
    else:
        # 赤道坐标系 - 使用改进的方法
        c1_min = c1_min_orig % 360
        c1_max = c1_max_orig % 360
        
        if lon_span >= 360.0:
            # 全经度范围，但赤纬有限制
            # 需要特殊处理，因为赤道坐标的"带状"在银道坐标中不是简单的带状
            pix_ids_in_range = get_pixels_in_equatorial_box(nside, 0, 360, c2_min, c2_max)
        else:
            # 使用改进的赤道坐标查询方法
            pix_ids_in_range = get_pixels_in_equatorial_box(nside, c1_min, c1_max, c2_min, c2_max)
    
    if len(pix_ids_in_range) == 0:
        raise ValueError('No data points in the selected area.')
    
    # 数据提取与计算
    df_selected = df_global.iloc[pix_ids_in_range]
    
    pix_ids, dust_values = hpix_Delta(df_selected, params['d_min'], params['d_max'])
    
    npix = hp.nside2npix(nside)
    dust_map_galactic = np.full(npix, hp.UNSEEN)
    dust_map_galactic[pix_ids] = dust_values
    
    # 平滑处理
    if params['smoothing_sigma'] > 0:
        dust_map_galactic = hp.smoothing(dust_map_galactic, sigma=np.radians(params['smoothing_sigma']))
    
    DEG_PER_PIXEL = 0.006; DPI = 300; MIN_PIXEL_DIM = 512
    MAX_PIXEL_DIM = 5120; AXES_TO_FIG_RATIO = 0.80; MAX_ASPECT_RATIO = 10.0
    
    if lon_span >= 360.0:
        range_lon = 360.0
    elif c1_min > c1_max:
        range_lon = (360.0 - c1_min) + c1_max
    else:
        range_lon = c1_max - c1_min
    
    range_lat = c2_max - c2_min
    if range_lon < 1e-6: range_lon = 1.0
    if range_lat < 1e-6: range_lat = 1.0
    
    ax_pixels_lon = range_lon / DEG_PER_PIXEL
    ax_pixels_lat = range_lat / DEG_PER_PIXEL
    aspect_ratio = ax_pixels_lon / ax_pixels_lat
    if aspect_ratio > MAX_ASPECT_RATIO: ax_pixels_lat = ax_pixels_lon / MAX_ASPECT_RATIO
    elif aspect_ratio < 1 / MAX_ASPECT_RATIO: ax_pixels_lon = ax_pixels_lat / MAX_ASPECT_RATIO
    if min(ax_pixels_lon, ax_pixels_lat) < MIN_PIXEL_DIM:
        scale = MIN_PIXEL_DIM / min(ax_pixels_lon, ax_pixels_lat)
        ax_pixels_lon *= scale; ax_pixels_lat *= scale
    if max(ax_pixels_lon, ax_pixels_lat) > MAX_PIXEL_DIM:
        scale = MAX_PIXEL_DIM / max(ax_pixels_lon, ax_pixels_lat)
        ax_pixels_lon *= scale; ax_pixels_lat *= scale
    fig_width  = (ax_pixels_lon / AXES_TO_FIG_RATIO) / DPI
    fig_height = (ax_pixels_lat / AXES_TO_FIG_RATIO) / DPI
    base_font_size = min(fig_height, fig_width) * 1.5
    dynamic_font_size = np.clip(base_font_size, 8, 24)
    
    # 绘图
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    if coord_system == 'equatorial':
        rot = hp.Rotator(coord=['G', 'C'], inv=False)
        map_to_plot = rot.rotate_map_pixel(dust_map_galactic)
        xlabel, ylabel = 'RA (°)', 'Dec (°)'
        lonra, latra = [c1_min_orig, c1_max_orig], [c2_min_orig, c2_max_orig]
    else:
        map_to_plot = dust_map_galactic
        xlabel, ylabel = 'Galactic Longitude l (°)', 'Galactic Latitude b (°)'
        lonra, latra = [c1_min_orig, c1_max_orig], [c2_min_orig, c2_max_orig]
    
    title_text = f"Dust map {params['d_min']:.2f}-{params['d_max']:.2f} kpc ({coord_system.capitalize()})"
    
    # 使用cartview绘制
    hp.cartview(
        map_to_plot, fig=fig.number, norm=params['norm'], format='%g',
        min=params['vmin'], max=params['vmax'], cmap=params['colormap'],
        lonra=lonra, latra=latra,
        title="", notext=True, cbar=False,
    )
    
    ax = plt.gca()
    ax.set_title(title_text, size=dynamic_font_size * 1.2, pad=dynamic_font_size * 0.8)
    
    # 颜色条
    if params['norm'] != 'hist':
        if params['norm'] == 'log':
            norm = mcolors.LogNorm(vmin=max(params['vmin'], 1e-9), vmax=params['vmax'])
        else:
            norm = mcolors.Normalize(vmin=params['vmin'], vmax=params['vmax'])
        sm = plt.cm.ScalarMappable(cmap=params['colormap'], norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.04, fraction=0.04, shrink=0.8)
        cbar.set_label('mag/kpc', size=dynamic_font_size)
        cbar.ax.tick_params(labelsize=dynamic_font_size * 0.9)
    
    # 网格和标签
    overlay_ax = fig.add_axes(ax.get_position(), frameon=False)
    overlay_ax.set_xlim(lonra[1], lonra[0])
    overlay_ax.set_ylim(latra[0], latra[1])
    overlay_ax.grid(True, color='black', linestyle='--', alpha=0.3)
    overlay_ax.set_xlabel(xlabel, size=dynamic_font_size)
    overlay_ax.set_ylabel(ylabel, size=dynamic_font_size)
    overlay_ax.tick_params(labelsize=dynamic_font_size * 0.9)
    
    # 保存图像
    filename = f"dust_map_{uuid.uuid4().hex}.png"
    save_path = str(Path(output_dir) / filename)
    plt.savefig(save_path, dpi=DPI, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    
    return {
        'path': save_path,
        'url': f"{public_base_url.rstrip('/')}/files/{filename}",
        'filename': filename
    }


def plot_partial_dust_map_ic_stilts(params, df_global, output_dir, public_base_url, stilts_command="stilts", timeout=300):
    """
    根据中心点和半径参数生成尘埃密度图 (STILTS版本)。
    该函数处理数据查询、计算、生成FITS文件，并调用STILTS绘图。
    
    Args:
        params (dict): 包含绘图参数的字典，必须包含:
            'lon_center', 'lat_center', 'radius_deg', 'coord_system',
            'd_min', 'd_max', 以及其他可选绘图参数。
    """
    nside = 1024
    npix = hp.nside2npix(nside)
    
    # 从统一后的参数中解包
    query_coord_system = params['coord_system']
    lon_center, lat_center = params['lon_center'], params['lat_center']
    radius_deg = params['radius_deg']
    d_min, d_max = params['d_min'], params['d_max']
    
    user_vmin = params.get('vmin', None)
    user_vmax = params.get('vmax', None)
    user_norm = params.get('norm', 'linear')
    user_colormap = params.get('colormap', 'viridis')
    should_mark_region = params.get('mark_region', False)
    bubble_diameter = params.get('bubble_diameter')
    mark_color = params.get('mark_color', 'black')
    show_color_bar = params.get('show_color_bar', True)
 
    is_full_sky = radius_deg >= 90.0
    expanded_radius_deg = min(radius_deg * 1.1, 180.0)
    query_radius_deg = 90.0

    if is_full_sky:
        pix_ids_in_range = np.arange(npix, dtype=np.int64)
    else:
        if query_coord_system == 'equatorial':
            galactic_center = SkyCoord(
                ra=lon_center * u.deg,
                dec=lat_center * u.deg,
                frame='icrs'
            ).galactic
            query_lon = galactic_center.l.deg
            query_lat = galactic_center.b.deg
        else:
            query_lon = lon_center
            query_lat = lat_center

        center_vec = hp.ang2vec(query_lon, query_lat, lonlat=True)
        pix_ids_in_range = np.asarray(
            hp.query_disc(
                nside,
                center_vec,
                radius=np.radians(query_radius_deg),
                nest=False,
                inclusive=True
            ),
            dtype=np.int64
        )

    if len(pix_ids_in_range) == 0:
        raise ValueError('No data points in the selected area.')

    df_selected = df_global if is_full_sky else df_global.iloc[pix_ids_in_range]
    pix_ids, dust_values = hpix_Delta(df_selected, d_min, d_max)
 
    dust_map_galactic = np.full(npix, hp.UNSEEN)
    if len(pix_ids) > 0:
        dust_map_galactic[pix_ids] = dust_values
 
    if params.get('smoothing_sigma', 0) > 0:
        sigma_rad = np.radians(params['smoothing_sigma'])
        valid_mask = (dust_map_galactic != hp.UNSEEN)
        dust_map_galactic[valid_mask] = hp.smoothing(dust_map_galactic, sigma=sigma_rad, lmax=2*nside)[valid_mask]
 
    valid_values = dust_map_galactic[dust_map_galactic != hp.UNSEEN]
    if len(valid_values) == 0:
        raise ValueError('No valid dust data in the selected distance range.')
    
    dust_map_galactic[dust_map_galactic == hp.UNSEEN] = np.nan
 
    with tempfile.TemporaryDirectory() as temp_dir:
        fits_input_path = os.path.join(temp_dir, 'dust_map_galactic_data.fits')
        png_output_path = os.path.join(temp_dir, 'dust_map_render.png')
 
        # --- FITS 文件生成 ---
        dust_map_nested = hp.reorder(dust_map_galactic, r2n=True)
        fits_col = fits.Column(name='I', format='D', array=dust_map_nested)
        hdu = fits.BinTableHDU.from_columns([fits_col])
        hdu.header['COORDSYS'] = ('G', 'Coordinate system: Galactic')
        hdu.header['PIXTYPE']  = 'HEALPix'
        hdu.header['ORDERING'] = 'NESTED'
        hdu.header['NSIDE']    = nside
        hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
        hdul.writeto(fits_input_path, overwrite=True)
 
        # --- STILTS 命令构建 ---
        plot_view_system = 'equatorial' if query_coord_system == 'equatorial' else 'galactic'
        
        stilts_colormap_map = {
            'viridis': 'viridis', 'plasma': 'plasma', 'inferno': 'inferno',
            'magma': 'magma', 'cividis': 'cividis',
            # 为 Spectral_r 提供一个合理的近似
            # 'spectral_r': '#9e0142-#d53e4f-#fdae61-#ffffbf-#66c2a5-#3288bd-#',
            'spectral_r': '#5e4fa2-#3288bd-#66c2a5-#ffffbf-#fdae61-#d53e4f-#9e0142',
            'pride': 'pride', 'sron': 'sron', 'rainbow': 'rainbow', 'iceburn': 'iceburn',
            'greys': 'greyscale'
        }
        stilts_colormap = stilts_colormap_map.get(user_colormap.lower(), 'viridis')
 
        stilts_norm_map = {'log': 'log', 'hist': 'histogram', 'linear': 'linear'}
        if user_norm.lower() in ['histogram', 'histolog', 'sqrt', 'square', 'acos', 'cos']:
            stilts_norm = user_norm.lower()
        else:
            stilts_norm = stilts_norm_map.get(user_norm.lower(), 'histogram')
 
        stilts_cmd_base = [
            stilts_command, '-batch', 'plot2sky', f'in={fits_input_path}', 
            f'title={d_min:.2f}-{d_max:.2f}[kpc]',
            'texttype=latex',
            'layer=healpix', 
            'datasys=galactic', 
            f'viewsys={plot_view_system}',
            'value=I', 
            
            'grid=true', 
            'labelpos=ExternalSys',
            'gridcolor=black',
            'labelcolor=black',
            # 'gridaa=true',      # 开启抗锯齿
            
            'sex=false', 
            'gridtrans=0',
            'legend=false',
        ]

        if show_color_bar:
            stilts_cmd_base.extend([
                'auxvisible=true',
                f'auxlabel=Density [mag/kpc]',
                f'auxmap={stilts_colormap}',
                f'auxfunc={stilts_norm}',
                'auxwidth=60',
                'auxcrowd=0.3',
            ])
        else:
            stilts_cmd_base.extend([
                'auxvisible=false',
                # f'auxlabel=Density [mag/kpc]',
                f'auxmap={stilts_colormap}',
                f'auxfunc={stilts_norm}',
                # 'auxwidth=60',
                # 'auxcrowd=0.3',
            ])

        # --- 动态设置字体和尺寸 ---
        if is_full_sky and not bubble_diameter:
            # 全天图的出版级字体和尺寸设置
            font_config = ['fontsize=72'] 
            image_config = [
                'projection=aitoff',
                # f'clon={lon_center}',
                # f'clat={lat_center}',
                'xpix=4000',
                'ypix=2000',
                # 'radius=180',
            ]
            stilts_cmd_extra = font_config + image_config
        else:
            # 局部天区图的出版级字体和尺寸设置
            font_config = ['fontsize=72'] 
            view_radius = expanded_radius_deg 
            image_config = [
                'projection=sin', f'clon={lon_center}', f'clat={lat_center}',
                f'radius={view_radius}', 'ypix=2000'
            ]
            if show_color_bar:
                image_config.append('xpix=2000')
            else:
                image_config.append('xpix=1700')
            stilts_cmd_extra = font_config + image_config
            
        stilts_cmd_base.extend(stilts_cmd_extra)
 
        # --- Vmin/Vmax 设置 ---
        if user_vmin is not None and user_vmax is not None and user_vmax > user_vmin:
            vmin, vmax = user_vmin, user_vmax
        else:
            vmin, vmax = np.nanmin(valid_values), np.nanmax(valid_values)
            if vmax - vmin < 1e-9: vmax = vmin + 1e-9
        stilts_cmd_base.extend([f'auxmin={vmin}', f'auxmax={vmax}'])
        
        stilts_cmd_extra = []
        # if is_full_sky:
        #     stilts_cmd_extra = ['projection=aitoff', 'xpix=3600', 'ypix=1800', 'fontSize=40']
        # else:
        #     # 给STILTS的半径稍微增加一点边距 (e.g., 10%)，避免图像边缘被裁切
        #     view_radius = radius_deg * 1.1 
        #     stilts_cmd_extra = [
        #         'projection=sin', f'clon={lon_center}', f'clat={lat_center}',
        #         f'radius={view_radius}', 'xpix=2000', 'ypix=2000', 'fontSize=40'
        #     ]
 
        if should_mark_region:
            if bubble_diameter is not None:
                mark_radius_deg = bubble_diameter / 2.0
            else:
                mark_radius_deg = radius_deg / 2.5

            area_frame = 'ICRS' if plot_view_system == 'equatorial' else 'GALACTIC'
            area_expr = f"\"Circle {area_frame} {lon_center} {lat_center} {mark_radius_deg}\""

            circle_csv_path = os.path.join(temp_dir, 'circle_data.csv')
            with open(circle_csv_path, 'w') as f:
                f.write("dummy\n")
                f.write("1\n")
            mark_layer_cmd = [
                'layer2=area',
                f'in2={circle_csv_path}',
                f'area2={area_expr}',
                'areatype2=STC-S',
                'polymode2=outline',
                f'color2={mark_color}',
                'thick2=3',
                'shading2=flat',
            ]
            stilts_cmd_extra.extend(mark_layer_cmd)

        stilts_cmd_final = stilts_cmd_base + stilts_cmd_extra + [f'out={png_output_path}']
 
        subprocess.run(stilts_cmd_final, text=True, check=True,
                       capture_output=True, cwd=temp_dir, timeout=timeout)

        unique_id = uuid.uuid4().hex[:8]
        filename = f"dust_map_{int(time.time())}_{unique_id}.png"
        
        save_dir = str(output_dir)
        os.makedirs(save_dir, exist_ok=True)
        final_save_path = os.path.join(save_dir, filename)
        
        shutil.move(png_output_path, final_save_path)
        
    return {
        'path': final_save_path,
        'url': f"{public_base_url.rstrip('/')}/files/{filename}",
        'filename': filename
    }

def plot_dust_slice_web(form_data, query_dust):
    """
    根据 Web 表单数据，在笛卡尔坐标系中绘制尘埃分布的切片或厚板，并返回图像的URL。
    此函数是 galactic_dust_plotter.py 的 Web 适配版本。
    """
    # --- 全局绘图样式设置 ---
    plt.rcParams.update({
        "font.size": 15, "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.titlesize": 18, "axes.labelsize": 16,
        "xtick.labelsize": 14, "ytick.labelsize": 14,
        "axes.linewidth": 1.2, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "axes.grid": True,
        "grid.alpha": 0.5, "grid.linestyle": "--", "grid.color": "gray",
    })
    # 1. 从 form_data 解包参数
    axis1, axis2, fixed_axis = form_data['axis1'].upper(), form_data['axis2'].upper(), form_data['fixed_axis'].upper()
    range1 = [form_data['range1_min'], form_data['range1_max']]
    range2 = [form_data['range2_min'], form_data['range2_max']]
    fixed_range = [form_data['fixed_range_min'], form_data['fixed_range_max']]
    
    resolution_pc = form_data['resolution_pc']
    smooth_sigma_pc = form_data['smooth_sigma_pc']
    vmin, vmax = form_data['vmin'], form_data['vmax']
    cmap = form_data['colormap']
    norm_str = form_data['norm']
    aggregate = form_data['aggregate']
    # 2. 单位转换与坐标轴映射
    lower_map = {'X': 'x', 'Y': 'y', 'Z': 'z'}
    ax1_lower, ax2_lower, fix_lower = lower_map[axis1], lower_map[axis2], lower_map[fixed_axis]
    resolution_kpc = resolution_pc / 1000.0
    smooth_sigma_kpc = smooth_sigma_pc / 1000.0
    # 3. 构建采样网格
    vals1 = np.arange(range1[0] + resolution_kpc/2, range1[1], resolution_kpc)
    vals2 = np.arange(range2[0] + resolution_kpc/2, range2[1], resolution_kpc)
    A, B = np.meshgrid(vals1, vals2, indexing='ij')
    n1, n2 = A.shape
    start, end = fixed_range
    if start == end:
        fixed_vals = np.array([start])
        title_slice = f"{fixed_axis} = {start:.3g} kpc"
    else:
        fixed_vals = np.arange(start + resolution_kpc/2, end, resolution_kpc)
        if len(fixed_vals) == 0: fixed_vals = np.array([(start + end) / 2])
        title_slice = f"{fixed_axis} in [{start:.3g}, {end:.3g}] kpc"
    n3 = len(fixed_vals)
    # 4. 创建三维坐标点云
    coords = {'x': None, 'y': None, 'z': None}
    coords[ax1_lower] = np.repeat(A[:, :, np.newaxis], n3, axis=2)
    coords[ax2_lower] = np.repeat(B[:, :, np.newaxis], n3, axis=2)
    fv_broadcast = fixed_vals[np.newaxis, np.newaxis, :]
    coords[fix_lower] = np.broadcast_to(fv_broadcast, (n1, n2, n3))
    X, Y, Z = coords['x'].ravel(), coords['y'].ravel(), coords['z'].ravel()
    # 5. 坐标转换并查询尘埃
    cart = CartesianRepresentation(X * u.kpc, Y * u.kpc, Z * u.kpc)
    skycoord = SkyCoord(cart, frame='galactic')
    l, b, d = skycoord.l.deg, skycoord.b.deg, skycoord.distance.kpc
    _, dust, _, max_d = query_dust(l, b, d)
    dust = np.asarray(dust, dtype=float)
    dust[d > np.asarray(max_d, dtype=float)] = np.nan
    # 6. 数据聚合与平滑
    dust_cube = dust.reshape(n1, n2, n3)
    agg_func = np.nanmean if aggregate == 'mean' else np.nanmedian
    agg2d = agg_func(dust_cube, axis=2)
    # 如果聚合后所有值都是NaN，说明区域内无数据
    if np.all(np.isnan(agg2d)):
        return {'url': None, 'error': 'No data found in the selected volume.'}
    # 高斯平滑（安全处理NaN）
    agg2d_filled = np.nan_to_num(agg2d, nan=0.0) # 用0填充NaN
    smooth_sigma_pixels = smooth_sigma_kpc / resolution_kpc
    smoothed = gaussian_filter(agg2d_filled, sigma=smooth_sigma_pixels)
    smoothed[np.isnan(agg2d)] = np.nan # 将原始NaN位置恢复为NaN
    # 7. 绘图
    Zplot = smoothed.T
    fig, ax = plt.subplots(figsize=(10, 8.5), constrained_layout=True)
    if norm_str == 'log':
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    
    im = ax.imshow(
        Zplot, extent=[vals1.min(), vals1.max(), vals2.min(), vals2.max()],
        origin='lower', aspect='equal', cmap=cmap, norm=norm
    )
    # 设置标签和标题
    ax.set_xlabel(f"{axis1} [kpc]")
    ax.set_ylabel(f"{axis2} [kpc]")
    title_agg = f"({aggregate})" if n3 > 1 else ""
    ax.set_title(f"Dust Density in {axis1}{axis2} Plane\n at {title_slice} {title_agg}")
    # 设置色条
    cbar = fig.colorbar(im, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("Dust Density [mag/kpc]", size=15)
    if isinstance(norm, LogNorm):
        ticks = np.geomspace(vmin, vmax, num=7)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{v:.2g}" for v in ticks])
    ax.set_aspect('equal', adjustable='box')
    # 8. 将图像保存到内存并编码为base64
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode('utf-8')
    plot_url = f'data:image/png;base64,{img_str}'

    return {'url': plot_url}


def load_dpr_superbubbles(superbubble_path):
    if not superbubble_path:
        return pd.DataFrame()
    table_path = Path(superbubble_path)
    if not table_path.is_file():
        raise FileNotFoundError("Superbubble table is unavailable")
    needed = {
        "id", "mark", "shape",
        "center_x_kpc", "center_y_kpc", "center_z_kpc",
        "a_radius_kpc", "b_radius_kpc", "c_radius_kpc", "angle_deg",
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_z_kpc", "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
    }
    df = pd.read_csv(table_path, usecols=lambda c: c in needed)
    for col in df.columns:
        if col != "shape":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def dpr_xy_terms(center_xy, axis_a, axis_b, angle_deg, unit_t, point_t0):
    if not np.isfinite(axis_a) or not np.isfinite(axis_b) or axis_a <= 0 or axis_b <= 0:
        return None
    phi = np.radians(float(angle_deg) if np.isfinite(angle_deg) else 0.0)
    major = np.array([np.cos(phi), np.sin(phi)])
    minor = np.array([-np.sin(phi), np.cos(phi)])
    delta = point_t0 - center_xy

    alpha_maj = float(np.dot(unit_t, major))
    alpha_min = float(np.dot(unit_t, minor))
    beta_maj = float(np.dot(delta, major))
    beta_min = float(np.dot(delta, minor))

    qa = 1.0 / (axis_a * axis_a)
    qb = 1.0 / (axis_b * axis_b)
    a_term = alpha_maj ** 2 * qa + alpha_min ** 2 * qb
    b_term = alpha_maj * beta_maj * qa + alpha_min * beta_min * qb
    c_term = beta_maj ** 2 * qa + beta_min ** 2 * qb
    if a_term <= DPR_EPS:
        return None
    return a_term, b_term, c_term


def dpr_xy_min_from_terms(terms):
    a_term, b_term, c_term = terms
    return c_term - b_term * b_term / a_term


def dpr_choose_slice_position(center_xy, axis_a, axis_b, angle, unit_t, unit_n, point_t0, half_width):
    def terms_at(offset):
        shifted_point_t0 = point_t0 + offset * unit_n
        terms = dpr_xy_terms(center_xy, axis_a, axis_b, angle, unit_t, shifted_point_t0)
        if terms is None:
            return None, np.inf
        return terms, dpr_xy_min_from_terms(terms)

    central_terms, central_xy_min = terms_at(0.0)
    if central_terms is not None and central_xy_min <= 1.0 + 1e-9:
        return central_terms, 0.0, central_xy_min
    if half_width <= DPR_EPS:
        return central_terms, 0.0, central_xy_min

    candidates = []
    for offset in (-half_width, half_width):
        terms, xy_min = terms_at(offset)
        candidates.append((xy_min, offset, terms))

    lo, hi = -half_width, half_width
    for _ in range(64):
        left = lo + (hi - lo) / 3.0
        right = hi - (hi - lo) / 3.0
        _, left_xy_min = terms_at(left)
        _, right_xy_min = terms_at(right)
        if left_xy_min < right_xy_min:
            hi = right
        else:
            lo = left
    offset = 0.5 * (lo + hi)
    terms, xy_min = terms_at(offset)
    candidates.append((xy_min, offset, terms))

    best_xy_min, best_offset, best_terms = min(candidates, key=lambda item: item[0])
    return best_terms, best_offset, best_xy_min


def dpr_compute_cut(row, unit_t, point_t0, unit_n, half_width=0.0):
    shape = str(row["shape"]).strip().lower()
    sb_id = int(row["id"])
    mark = int(row["mark"]) if np.isfinite(row["mark"]) else -1
    cap_section_ids = DPR_SB21_CAP_SECTION_IDS | DPR_DEFAULT_EXTRA_CAP_SECTION_IDS

    if shape == "cylinder":
        center_xy = np.array([row["center_x_kpc"], row["center_y_kpc"]], dtype=float)
        axis_a = float(row["a_radius_kpc"])
        axis_b = float(row["b_radius_kpc"])
        angle = float(row["angle_deg"]) if np.isfinite(row["angle_deg"]) else 0.0
        z_center = float(row["center_z_kpc"])
        z_radius = DPR_CYLINDER_HALF_HEIGHT_KPC
        cut_kind = "rectangle"
        section_model = "cylinder"
    elif sb_id in cap_section_ids:
        center_xy = np.array([
            row.get("xy_plane_center_x_kpc", row["center_x_kpc"]),
            row.get("xy_plane_center_y_kpc", row["center_y_kpc"]),
        ], dtype=float)
        axis_a = float(row["xy_plane_a_kpc"])
        axis_b = float(row["xy_plane_b_kpc"])
        angle = float(row["xy_plane_angle_deg"]) if np.isfinite(row["xy_plane_angle_deg"]) else float(row["angle_deg"])
        xy_z = float(row["xy_plane_z_kpc"])
        c_full = float(row["c_radius_kpc"])
        z_apex = float(row["center_z_kpc"]) - c_full if mark == 1 else float(row["center_z_kpc"]) + c_full
        z_center = xy_z
        z_radius = abs(xy_z - z_apex)
        cut_kind = "ellipse_arc"
        section_model = "ellipsoid_cap"
    else:
        center_xy = np.array([row["center_x_kpc"], row["center_y_kpc"]], dtype=float)
        axis_a = float(row["a_radius_kpc"])
        axis_b = float(row["b_radius_kpc"])
        angle = float(row["angle_deg"]) if np.isfinite(row["angle_deg"]) else 0.0
        z_center = float(row["center_z_kpc"])
        z_radius = float(row["c_radius_kpc"])
        cut_kind = "ellipse_arc"
        section_model = "full_ellipsoid"

    terms, slice_offset, xy_min = dpr_choose_slice_position(
        center_xy, axis_a, axis_b, angle, unit_t, unit_n, point_t0, half_width
    )
    if terms is None:
        return None
    a_term = terms[0]
    t_center = -terms[1] / a_term
    if xy_min > 1.0 + 1e-9:
        return None
    remaining = max(0.0, 1.0 - xy_min)
    t_radius = float(np.sqrt(remaining / a_term))
    if cut_kind == "ellipse_arc":
        z_radius = z_radius * np.sqrt(remaining)
    if t_radius <= DPR_EPS or z_radius <= DPR_EPS:
        return None

    result = {
        "id": sb_id,
        "mark": mark,
        "cut_kind": cut_kind,
        "section_model": section_model,
        "t_center_kpc": t_center,
        "z_center_kpc": z_center,
        "t_radius_kpc": t_radius,
        "z_radius_kpc": z_radius,
    }
    if section_model == "full_ellipsoid":
        result.update({
            "full_t_center_kpc": t_center,
            "full_z_center_kpc": z_center,
            "full_t_radius_kpc": t_radius,
            "full_z_radius_kpc": z_radius,
        })
    return result


def dpr_compute_all_cuts(superbubbles, unit_t, point_t0, unit_n):
    rows = []
    for _, row in superbubbles.iterrows():
        cut = dpr_compute_cut(row, unit_t, point_t0, unit_n)
        if cut is not None:
            rows.append(cut)
    return pd.DataFrame(rows)


def dpr_ellipsoid_section_params(cut):
    if cut["cut_kind"] == "rectangle":
        return (
            cut["t_center_kpc"], cut["z_center_kpc"],
            cut["t_radius_kpc"], cut["z_radius_kpc"],
        )
    full_values = [
        cut.get("full_t_center_kpc", np.nan),
        cut.get("full_z_center_kpc", np.nan),
        cut.get("full_t_radius_kpc", np.nan),
        cut.get("full_z_radius_kpc", np.nan),
    ]
    if all(np.isfinite(v) for v in full_values):
        return tuple(full_values)
    return (
        cut["t_center_kpc"], cut["z_center_kpc"],
        cut["t_radius_kpc"], cut["z_radius_kpc"],
    )


def draw_dpr_cut(ax, cut, color="black"):
    tc, zc, tr, zr = dpr_ellipsoid_section_params(cut)
    if cut["cut_kind"] == "rectangle":
        ax.add_patch(Rectangle(
            (tc - tr, zc - 2 * zr),
            width=2 * tr, height=4 * zr,
            edgecolor=color, facecolor="none", linewidth=1.0, alpha=0.9, zorder=6,
        ))
        return
    mark = int(cut["mark"])
    if mark == 1:
        angles = np.linspace(np.pi, 2 * np.pi, 180)
    elif mark == 2:
        angles = np.linspace(0, np.pi, 180)
    else:
        angles = np.linspace(0, 2 * np.pi, 240)
    t = tc + tr * np.cos(angles)
    z = zc + zr * np.sin(angles)
    ax.plot(t, z, color=color, linewidth=1.1, alpha=0.95, zorder=6)


def draw_dpr_cuts(ax, cuts):
    if cuts.empty:
        return
    for _, cut in cuts.iterrows():
        color = "black"
        draw_dpr_cut(ax, cut, color)
        label_t, label_z, _, label_z_radius = dpr_ellipsoid_section_params(cut)
        mark = int(cut["mark"])
        if mark == 1:
            label_z -= 0.35 * label_z_radius
        elif mark == 2:
            label_z += 0.35 * label_z_radius
        ax.text(
            label_t, label_z, f"SB{int(cut['id'])}",
            fontsize=DPR_SB_LABEL_FONTSIZE, fontweight="bold", ha="center", va="center",
            color=color, zorder=7,
            path_effects=[pe.withStroke(linewidth=1.4, foreground="white")],
        )


def plot_oblique_slice_web(form_data, df3d, superbubble_path):
    """
    只生成切片内尘埃垂直分布图（右图）
    特性：
    1. 物理分辨率固定为 0.015 kpc/pixel
    2. 像素为正方形 (Square Pixels)
    3. T/Z Range 仅起到"相机变焦"作用，不改变底层数据网格
    """
    # 全局绘图样式设置
    plt.rcParams.update({
        "font.size": 12,
        "font.family": "DejaVu Sans",
        "axes.unicode_minus": False,
        "axes.linewidth": 1.2,
    })
    
    
    # 解包参数
    angle_degrees = form_data['angle_degrees']
    offset = form_data['offset']
    half_width = form_data['half_width']
    theta = np.radians(angle_degrees)
    

    dir_vec = np.array([np.cos(theta), np.sin(theta)])
    normal_vec = np.array([-np.sin(theta), np.cos(theta)])
    

    x_values = df3d['X'].to_numpy()
    y_values = df3d['Y'].to_numpy()
    z_values = df3d['Z'].to_numpy()
    dust_values = df3d['dust'].to_numpy()

    distances = (x_values * normal_vec[0] +
                 y_values * normal_vec[1] -
                 offset)

    # 筛选切片带内的数据
    mask = np.abs(distances) <= half_width

    if not np.any(mask):
        return {'url_right': None, 'error': '切片带内无数据'}


    x0 = offset * normal_vec[0]
    y0 = offset * normal_vec[1]
    point_t0 = np.array([x0, y0])

    # 准备s, z, dust数据
    s = ((x_values[mask] - point_t0[0]) * dir_vec[0] +
         (y_values[mask] - point_t0[1]) * dir_vec[1])
    z = z_values[mask]
    dust = dust_values[mask]
    
    # 移除NaN值
    mask_valid = ~np.isnan(dust)
    s, z, dust = s[mask_valid], z[mask_valid], dust[mask_valid]
    
    if len(s) == 0:
        return {'url_right': None, 'error': '切片带内无有效尘埃数据'}
    
    # 设定固定的空间分辨率 (0.015 kpc/pixel)
    PIXEL_RES = 0.02
    
    # 设定全局物理边界 (Hardcode)
    # 这个范围应足够大，以覆盖所有可能的数据区域
    # 即使数据只有 -2到2，我们也生成 -5到5的网格，这保证了网格对齐的一致性
    GLOBAL_S_MIN, GLOBAL_S_MAX = -5.0, 5.0
    GLOBAL_Z_MIN, GLOBAL_Z_MAX = -3.0, 3.0
    
    # 自动生成边缘数组
    # arange 保证步长严格一致
    s_edges = np.arange(GLOBAL_S_MIN, GLOBAL_S_MAX + PIXEL_RES, PIXEL_RES)
    z_edges = np.arange(GLOBAL_Z_MIN, GLOBAL_Z_MAX + PIXEL_RES, PIXEL_RES)
    
    hist, _, _ = np.histogram2d(s, z, bins=[s_edges, z_edges], weights=dust)
    counts, _, _ = np.histogram2d(s, z, bins=[s_edges, z_edges])
    mean_dust = np.divide(hist, counts, out=np.full_like(hist, np.nan), where=counts != 0).T
    
    # 获取用户设置的显示范围 (Zoom window)
    display_s_range = [form_data['s_range_min'], form_data['s_range_max']]
    display_z_range = [form_data['z_range_min'], form_data['z_range_max']]
    
    # 计算可视区域的物理尺寸
    view_width = display_s_range[1] - display_s_range[0]
    view_height = display_z_range[1] - display_z_range[0]
    
    # 物理长宽比
    view_aspect_ratio = view_width / view_height
    
    # 限制画布长宽比 (防止极端长条图)
    limited_aspect_ratio = np.clip(view_aspect_ratio, 0.25, 4.0)
    
    # 计算画布像素尺寸 (保证最短边至少 6 英寸)
    base_size = 6.0
    if limited_aspect_ratio >= 1.0:
        # 横向图
        fig_height = base_size
        fig_width = base_size * limited_aspect_ratio
    else:
        # 纵向图
        fig_width = base_size
        fig_height = base_size / limited_aspect_ratio
        
    # 确保最小尺寸兜底
    fig_width = max(fig_width, 6.0)
    fig_height = max(fig_height, 6.0)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # 裁剪数据范围（避免对数scale的警告）
    mean_dust_clipped = np.clip(mean_dust, form_data['vmin'], form_data['vmax'])
    
    mesh = ax.pcolormesh(
        s_edges, z_edges, mean_dust_clipped,
        cmap=form_data['colormap'],
        shading='flat',
        norm=LogNorm(vmin=form_data['vmin'], vmax=form_data['vmax'])
    )

    superbubbles = load_dpr_superbubbles(superbubble_path)
    cuts = dpr_compute_all_cuts(superbubbles, dir_vec, point_t0, normal_vec) if not superbubbles.empty else pd.DataFrame()
    draw_dpr_cuts(ax, cuts)

    ax.set_aspect('equal', adjustable='box')
    
    # 这步操作相当于相机的变焦，不会改变像素的物理大小
    ax.set_xlim(display_s_range)
    ax.set_ylim(display_z_range)
    
    # 标签与装饰
    ax.set_xlabel('T [kpc]', fontsize=13)
    ax.set_ylabel('Z [kpc]', fontsize=13)
    
    if abs(np.cos(theta)) < 1e-5:
        intercept_str = "b → ∞"
    else:
        intercept_val_kpc = offset / np.cos(theta)
        intercept_str = f"b = {round(1000 * intercept_val_kpc, 1):g} pc"

    ax.set_title(
        f'θ = {round(form_data["angle_degrees"], 1):g}°; {intercept_str}; Offset = {round(1000*form_data["offset"], 1):g} pc',
        fontsize=14, pad=12
    )
    ax.grid(True, linestyle='--', alpha=0.4)
    
    # 自适应颜色条高度
    # 如果图像很宽，缩小颜色条的aspect（让它更矮）
    if limited_aspect_ratio > 2.0:
        cbar_shrink = 0.7
        cbar_aspect = 20
    elif limited_aspect_ratio < 0.5:
        cbar_shrink = 0.95
        cbar_aspect = 35
    else:
        cbar_shrink = 0.85
        cbar_aspect = 25
    
    cbar = fig.colorbar(mesh, ax=ax, shrink=cbar_shrink, aspect=cbar_aspect)
    cbar.set_label('Dust [mag/kpc]', fontsize=12)
    
    if form_data.get('show_markers', False):
        ax.scatter([0], [0], color='red', marker='+', s=100, label='Origin', zorder=10)

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode('utf-8')
    
    return {'url_right': f'data:image/png;base64,{img_str}'}

def annulus_pixels(
    l_deg: float,
    b_deg: float,
    r1_deg: float = 1.0,
    r2_deg: float = 2.0,
    nside: int = 1024,
    nest: bool = False,
):
    """
    返回 nside 下中心在 (l, b) 的 r1~r2 度环形区域内的像素索引。
    l_deg: 银河经度 (deg)
    b_deg: 银河纬度 (deg)
    r1_deg, r2_deg: 内/外半径 (deg)，要求 r2 > r1
    """
    if r2_deg <= r1_deg:
        raise ValueError("r2 must be larger than r1")

    # healpy 使用球坐标：theta=90°-b（余纬度），phi=l，经度
    theta = np.deg2rad(90.0 - b_deg)
    phi = np.deg2rad(l_deg)
    vec = hp.ang2vec(theta, phi)

    r1_rad = np.deg2rad(r1_deg)
    r2_rad = np.deg2rad(r2_deg)

    # 外圆与内圆；均 inclusive，可避免空洞
    outer_pix = set(hp.query_disc(nside, vec, r2_rad, inclusive=True, nest=nest))
    inner_pix = set(hp.query_disc(nside, vec, r1_rad, inclusive=True, nest=nest))
    annulus_pix = outer_pix - inner_pix
    return sorted(annulus_pix)

from typing import Optional

def dust_profile_mean(
    l_center_deg: float,
    b_center_deg: float,
    r1_deg: float,
    r2_deg: float,
    nside: int = 1024,
    d_start: float = 2.0,
    d_end: float = 3.0,
    d_samples: int = 100,
    nest: bool = False,
    inner_radius_deg: Optional[float] = None,
    query_dust=None,
    max_points=2_000_000,
):
    """
    计算随距离 d 的 dust 平均值：外环 (r1~r2) 与内圆 (r<=inner_radius)。
    inner_radius_deg 为 None 时，默认使用 (2/3)*r1。
    返回距离数组、外环均值数组、内圆均值数组。
    """
    # 像素索引
    annulus_pix = annulus_pixels(l_center_deg, b_center_deg, r1_deg, r2_deg, nside, nest=nest)
    theta_inner = np.deg2rad(90.0 - b_center_deg)
    phi_inner = np.deg2rad(l_center_deg)
    vec_center = hp.ang2vec(theta_inner, phi_inner)
    inner_radius = inner_radius_deg if inner_radius_deg is not None else (2.0 / 3.0) * r1_deg
    inner_pix = hp.query_disc(nside, vec_center, np.deg2rad(inner_radius), inclusive=True, nest=nest)

    # 转为角度 (l, b)
    ann_theta, ann_phi = hp.pix2ang(nside, annulus_pix, nest=nest)
    inner_theta, inner_phi = hp.pix2ang(nside, inner_pix, nest=nest)
    ann_l = np.rad2deg(ann_phi)
    ann_b = 90.0 - np.rad2deg(ann_theta)
    inner_l = np.rad2deg(inner_phi)
    inner_b = 90.0 - np.rad2deg(inner_theta)

    d_vals = np.linspace(d_start, d_end, num=d_samples)

    # 向量化：预先铺开所有距离-像素组合，一次性查询 dustmaps3d
    ann_size = ann_l.size
    inner_size = inner_l.size

    if (ann_size + inner_size) * d_samples > max_points:
        raise ValueError("Bubble region exceeds the maximum sample count")
    l_ann = np.tile(ann_l, d_samples)
    b_ann = np.tile(ann_b, d_samples)
    d_ann = np.repeat(d_vals, ann_size)
    _, dust_ann_flat, _, _ = query_dust(l_ann, b_ann, d_ann)
    dust_ann = np.asarray(dust_ann_flat, dtype=float).reshape(d_samples, ann_size)
    ann_means = np.nanmean(dust_ann, axis=1)

    l_inn = np.tile(inner_l, d_samples)
    b_inn = np.tile(inner_b, d_samples)
    d_inn = np.repeat(d_vals, inner_size)
    _, dust_inn_flat, _, _ = query_dust(l_inn, b_inn, d_inn)
    dust_inn = np.asarray(dust_inn_flat, dtype=float).reshape(d_samples, inner_size)
    inner_means = np.nanmean(dust_inn, axis=1)

    deltas = ann_means - inner_means

    return d_vals, ann_means, inner_means, deltas


def half_width_main_peak(
    d_vals: np.ndarray,
    y_vals: np.ndarray,
    d_low: float,
    d_up: float,
    prominence_frac: float = 0.1,
):
    """寻找主峰的半高宽，返回峰位置、左右交点与宽度。"""
    mask = (d_vals >= d_low) & (d_vals <= d_up)
    x_win = d_vals[mask]
    y_win = y_vals[mask]
    if x_win.size == 0:
        raise ValueError("d_low~d_up 范围内没有采样点，无法标记峰值")

    y_min = np.nanmin(y_win)
    y_max = np.nanmax(y_win)
    if not np.isfinite(y_min) or not np.isfinite(y_max):
        return np.nan, np.nan, np.nan, np.nan, np.nan
    baseline = 0.0  # use zero as reference for FWHM
    amplitude = y_max - baseline
    if amplitude <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    prom = max(1e-9, amplitude * prominence_frac)
    distance = max(1, x_win.size // 20)
    peaks, props = find_peaks(y_win, prominence=prom, distance=distance)
    if peaks.size == 0:
        peak_idx = int(np.nanargmax(y_win))
    else:
        prominences = props.get("prominences", np.zeros_like(peaks, dtype=float))
        order = np.lexsort((-y_win[peaks], -prominences))
        peak_idx = int(peaks[order[0]])

    peak_x = float(x_win[peak_idx])
    peak_y = float(y_win[peak_idx])
    height = peak_y - baseline
    if not np.isfinite(height) or height <= 0:
        return peak_x, peak_y, np.nan, np.nan, np.nan
    half_y = baseline + 0.5 * height

    def interp_cross(start_idx: int, step: int, boundary: float):
        i = start_idx
        while 0 <= i + step < len(x_win):
            y0, y1 = y_win[i], y_win[i + step]
            x0, x1 = x_win[i], x_win[i + step]
            if not (np.isfinite(y0) and np.isfinite(y1)):
                i += step
                continue
            if (y0 - half_y) * (y1 - half_y) <= 0 and y0 != y1:
                t = (half_y - y0) / (y1 - y0)
                return x0 + t * (x1 - x0)
            i += step
        return boundary

    left = interp_cross(peak_idx, -1, x_win[0])
    right = interp_cross(peak_idx, 1, x_win[-1])
    width = right - left
    return peak_x, peak_y, left, right, width


def analyze_dust(
    l: float,
    b: float,
    d: float,
    d_low: float,
    d_up: float,
    diameter: float,
    inner_factor: float = 0.25,
    annulus_inner_factor: float = 0.375,
    annulus_outer_factor: float = 0.625,
    query_dust=None,
) -> str:
    """
    分析指定位置的尘埃分布，返回 base64 编码的图片。

    参数:
        l: 银河经度 (deg)
        b: 银河纬度 (deg)
        d: 参考距离 (kpc)，暂保留
        d_low: 距离范围下界 (kpc)
        d_up: 距离范围上界 (kpc)
        diameter: 直径 (deg)
        inner_factor: 内圆半径 = inner_factor * diameter
        annulus_inner_factor: 外环内径 = annulus_inner_factor * diameter
        annulus_outer_factor: 外环外径 = annulus_outer_factor * diameter

    返回:
        base64 编码的 PNG 图片字符串
    """
    # 计算区域尺寸
    inner_disk_radius_deg = inner_factor * diameter
    annulus_inner_deg = annulus_inner_factor * diameter
    annulus_outer_deg = annulus_outer_factor * diameter

    d_span = d_up - d_low
    d_start = d_low - d_span
    d_end = d_up + d_span

    # 计算 dust 随距离的外环与内圆平均值
    d_vals, ann_mean, inner_mean, delta_mean = dust_profile_mean(
        l,
        b,
        r1_deg=annulus_inner_deg,
        r2_deg=annulus_outer_deg,
        nside=1024,
        d_start=d_start,
        d_end=d_end,
        d_samples=100,
        inner_radius_deg=inner_disk_radius_deg,
        query_dust=query_dust,
    )

    delta_smooth = delta_mean

    # 高斯拟合
    def gaussian(x, amp, mu, sigma_g):
        return amp * np.exp(-0.5 * ((x - mu) / sigma_g) ** 2)

    fit_mask = (d_vals >= d_low) & (d_vals <= d_up) & np.isfinite(delta_mean)
    x_fit = d_vals[fit_mask]
    y_fit = delta_mean[fit_mask]
    half_level = np.nan
    err_low = err_high = np.nan

    gaussian_curve = None
    peak_d = peak_val = sigma = sigma_low = sigma_high = np.nan

    if x_fit.size >= 5:
        prom_guess = max(1e-9, (np.nanmax(y_fit) - np.nanmin(y_fit)) * 0.05)
        pk_idx, pk_props = find_peaks(y_fit, prominence=prom_guess)
        if pk_idx.size > 0:
            prominences = pk_props.get("prominences", np.zeros_like(pk_idx, dtype=float))
            best = int(pk_idx[np.nanargmax(prominences)])
        else:
            best = int(np.nanargmax(y_fit))

        mu0 = float(x_fit[best]) if np.isfinite(x_fit[best]) else float(np.nanmean(x_fit))
        span = max(0.05 * (d_up - d_low), (d_up - d_low) / 10.0)
        fit_local_mask = (x_fit >= mu0 - span) & (x_fit <= mu0 + span)
        x_use = x_fit[fit_local_mask] if np.count_nonzero(fit_local_mask) >= 5 else x_fit
        y_use = y_fit[fit_local_mask] if np.count_nonzero(fit_local_mask) >= 5 else y_fit

        amp0 = max(1e-6, np.nanmax(y_use))
        mu0 = x_use[np.nanargmax(y_use)] if np.isfinite(np.nanmax(y_use)) else float(np.nanmean(x_use))
        sigma0 = max(1e-3, span / 2.0)
        try:
            popt, pcov = curve_fit(
                gaussian,
                x_use,
                y_use,
                p0=[amp0, mu0, sigma0],
                bounds=([0, d_low, 1e-6], [np.inf, d_up, (d_up - d_low)]),
                maxfev=6000,
            )
            amp_f, mu_f, sigma_g = popt
            fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma_g
            peak_d = float(mu_f)
            peak_val = float(amp_f)
            sigma = float(fwhm)
            err_low = err_high = float(fwhm / 2.0)
            sigma_low = peak_d - err_low
            sigma_high = peak_d + err_high
            half_level = float(0.5 * amp_f)
            x_dense = np.linspace(d_vals.min(), d_vals.max(), 200)
            gaussian_curve = (x_dense, gaussian(x_dense, *popt))
        except (ValueError, RuntimeError):
            pass

    # 高斯拟合失败时回退半高宽方法
    if gaussian_curve is None:
        peak_d, peak_val, sigma_low, sigma_high, sigma = half_width_main_peak(
            d_vals, delta_smooth, d_low, d_up, prominence_frac=0.1
        )
        err_low = peak_d - sigma_low
        err_high = sigma_high - peak_d
        peak_idx_plot = int(np.nanargmin(np.abs(d_vals - peak_d))) if np.isfinite(peak_d) else None
        half_level = 0.5 * delta_mean[peak_idx_plot] if peak_idx_plot is not None and np.isfinite(delta_mean[peak_idx_plot]) else np.nan

    # 绘图
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    fig.set_facecolor("white")

    ax_left, ax_right = axes
    colors = {"ann": "#1f77b4", "inner": "#2ca02c", "delta": "#9467bd", "delta_raw": "#7f7f7f", "delta_smooth": "#17becf"}
    delta_style_left = dict(color=colors["delta"], linewidth=1.6, linestyle="--", marker="o", markersize=3, alpha=0.85)

    ax_left.plot(d_vals, ann_mean, label="Annulus dust mean", color=colors["ann"], linewidth=1.8, marker="o", markersize=3, alpha=0.9)
    ax_left.plot(d_vals, inner_mean, label="Inner disk dust mean", color=colors["inner"], linewidth=1.8, marker="o", markersize=3, alpha=0.9)
    ax_left.plot(d_vals, delta_mean, label="Δ (annulus - inner)", **delta_style_left)
    if np.isfinite(peak_d):
        ax_left.axvline(peak_d, color="#d62728", linestyle="--", linewidth=1.0, label="Fit peak d")
    ax_left.set_xlabel("Distance [kpc]")
    ax_left.set_ylabel("Dust [mag/kpc]")
    ax_left.set_title("Dust mean vs distance", pad=8)
    ax_left.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)
    for spine in ax_left.spines.values():
        spine.set_linewidth(0.8)
        spine.set_visible(True)
    ax_left.tick_params(labelsize=9)
    ax_left.legend(frameon=True, fontsize=9, framealpha=0.9, edgecolor="#cccccc")

    ax_right.plot(d_vals, delta_mean, label="Δ (raw)", alpha=0.6, linestyle="--", marker="o", markersize=3, color=colors["delta"], linewidth=1.0)
    if gaussian_curve is not None:
        x_dense, y_dense = gaussian_curve
        ax_right.plot(x_dense, y_dense, label="Gaussian fit", color="#1f78b4", linewidth=1.6, linestyle="-")
    if np.isfinite(err_low) and np.isfinite(err_high):
        peak_label = rf"Peak $d={peak_d:.2f}\pm{err_high:.2f}\,\mathrm{{kpc}}$"
    else:
        peak_label = f"Peak d={peak_d:.2f} kpc"
    ax_right.axvline(peak_d, color="#d62728", linestyle="--", linewidth=1.2, label=peak_label)
    # 添加 d_low 和 d_up 参考线
    ax_right.axvline(d_low, color="#2ca02c", linestyle="--", linewidth=1.0, label=f"d_low={d_low:.2f} kpc")
    ax_right.axvline(d_up, color="#2ca02c", linestyle="--", linewidth=1.0, label=f"d_up={d_up:.2f} kpc")
    if np.isfinite(peak_d) and np.isfinite(peak_val):
        ax_right.plot(peak_d, peak_val, marker="o", color="#d62728", markersize=5, label=f"Peak value {peak_val:.2f} mag/kpc")
    if np.isfinite(sigma):
        ax_right.axvspan(
            sigma_low,
            sigma_high,
            color="#d62728",
            alpha=0.12,
            label=f"sigma = {sigma:.2f} kpc",
        )
    if np.isfinite(sigma) and sigma > 0:
        half_span = sigma
        x_left = peak_d - half_span
        x_right = peak_d + half_span
        x_left = min(x_left, d_low)
        x_right = max(x_right, d_up)
        x_left = max(d_vals.min(), x_left)
        x_right = min(d_vals.max(), x_right)
        ax_right.set_xlim(x_left, x_right)
    else:
        ax_right.set_xlim(d_low, d_up)
    ax_right.set_xlabel("Distance [kpc]")
    ax_right.set_ylabel("Δ dust [mag/kpc]")
    ax_right.set_title("Peak localization & FWHM", pad=8)
    ax_right.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)
    ax_right.set_ylim(ax_left.get_ylim())

    # 添加 sigma 区间边界的文字标注
    if np.isfinite(sigma_low) and np.isfinite(sigma_high):
        x0, x1 = ax_right.get_xlim()
        y0, y1 = ax_right.get_ylim()
        x_pad = 0.01 * (x1 - x0)
        y_pad = 0.04 * (y1 - y0)
        span_y = y1 - y_pad
        ax_right.text(
            sigma_low + x_pad,
            span_y,
            f"{sigma_low:.2f}",
            color="#d62728",
            rotation=0,
            va="top",
            ha="left",
            fontweight="bold",
        )
        ax_right.text(
            sigma_high - x_pad,
            span_y,
            f"{sigma_high:.2f}",
            color="#d62728",
            rotation=0,
            va="top",
            ha="right",
            fontweight="bold",
        )

    for spine in ax_right.spines.values():
        spine.set_linewidth(0.8)
        spine.set_visible(True)
    ax_right.tick_params(labelsize=9)
    ax_right.legend(frameon=True, fontsize=9, framealpha=0.9, edgecolor="#cccccc")

    # 保存到内存并编码为 base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")

    return img_base64

def plot_schematic(
    diameter: float,
    inner_factor: float = 0.25,
    annulus_inner_factor: float = 0.375,
    annulus_outer_factor: float = 0.625,
) -> str:
    """
    绘制环形区域示意图，返回 base64 编码的图片。

    参数:
        diameter: 直径 (deg)
        inner_factor: 内圆半径 = inner_factor * diameter
        annulus_inner_factor: 外环内径 = annulus_inner_factor * diameter
        annulus_outer_factor: 外环外径 = annulus_outer_factor * diameter

    返回:
        base64 编码的 PNG 图片字符串
    """
    # 计算尺寸
    inner_r = inner_factor * diameter
    annulus_r_inner = annulus_inner_factor * diameter
    annulus_r_outer = annulus_outer_factor * diameter
    d_circle_r = diameter / 2

    fig, ax = plt.subplots(figsize=(6, 6))

    # 外环使用两个同心圆模拟：填充外圆，挖空内圆
    outer_ring = Circle((0, 0), annulus_r_outer, facecolor="#cfe2ff", edgecolor="none", alpha=0.85)
    inner_hole = Circle((0, 0), annulus_r_inner, facecolor="white", edgecolor="none")

    # 内圆
    inner_circle = Circle((0, 0), inner_r, facecolor="#9acd32", edgecolor="none")

    # 直径 D 的黑色圆
    d_circle = Circle((0, 0), d_circle_r, facecolor="none", edgecolor="black", lw=1.8, ls="-")

    ax.add_patch(outer_ring)
    ax.add_patch(inner_hole)
    ax.add_patch(inner_circle)
    ax.add_patch(d_circle)

    # 设定视图范围
    margin = diameter * 0.3
    limit = annulus_r_outer + margin
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", "box")

    # 去掉刻度与框线以突出示意
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()

    # 保存到内存并编码为 base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")

    return img_base64
