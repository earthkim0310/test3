# Streamlit RV (Radial Velocity) Exoplanet Simulator — Cloud friendly
# - No GUI backends (TkAgg 제거), matplotlib Agg 사용
# - Streamlit 위젯로 제어, st.pyplot으로 렌더링
# - 원형 궤도(e=0) 기본, 필요 시 e/ω 확장 가능
#
# 요구 패키지 (requirements.txt 예시)
# streamlit
# matplotlib
# numpy
# pandas
# pillow

import math
import numpy as np
import streamlit as st

import matplotlib
matplotlib.use("Agg")  # 서버/클라우드 환경
import matplotlib.pyplot as plt

# ----------------------------- Constants -----------------------------
G = 6.67430e-11  # m^3 kg^-1 s^-2
MSUN = 1.98847e30  # kg
MJUP = 1.89813e27  # kg
DAY = 86400.0  # s
C = 299_792_458.0  # m/s

# ----------------------------- Page Config -----------------------------
st.set_page_config(
    page_title="RV Exoplanet Simulator (Streamlit)",
    layout="wide",
)

st.title("🌌 Radial Velocity Exoplanet Simulator (Streamlit)")
st.caption(
    "Streamlit Cloud 호환 버전 — matplotlib 'Agg' 백엔드 사용, Tkinter/GUI 의존성 제거"
)

# ----------------------------- Sidebar Controls -----------------------------
with st.sidebar:
    st.header("⚙️ 파라미터")
    Ms = st.number_input("항성 질량 Ms (M☉)", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
    Mp = st.number_input("행성 질량 Mp (M♃)", min_value=0.01, max_value=20.0, value=1.0, step=0.01)
    P_days = st.number_input("공전 주기 P (일)", min_value=0.1, max_value=5000.0, value=365.25, step=0.1)
    inc_deg = st.slider("기울기 i (°)", min_value=0, max_value=90, value=90)
    phase = st.slider("위상 (0→1)", min_value=0.0, max_value=1.0, value=0.0, step=0.01)
    e = st.slider("이심률 e", min_value=0.0, max_value=0.9, value=0.0, step=0.01)
    omega_deg = st.slider("근일점 인수 ω (°) — e>0일 때 의미", min_value=0, max_value=360, value=90)

    st.markdown("---")
    t_window_factor = st.select_slider("시간창 길이 (주기의 배수)", options=[0.5, 1.0, 2.0, 3.0], value=1.0)
    npts = st.slider("샘플 수", min_value=200, max_value=4000, value=1000, step=100)

# ----------------------------- Derived -----------------------------
Ms_kg = Ms * MSUN
Mp_kg = Mp * MJUP
P = P_days * DAY
inc = math.radians(inc_deg)
omega = math.radians(omega_deg)

# 궤도반지름 a (케플러 3법칙; 이심률 무관): a^3 = G(Ms+Mp) P^2 / (4π^2)
a = (G * (Ms_kg + Mp_kg) * P * P / (4.0 * math.pi**2)) ** (1.0 / 3.0)  # meters

# ----------------------------- RV model -----------------------------
# 일반적인 RV 식 (표준 형식):
# v_r(t) = K [cos(θ(t) + ω) + e cos ω]  (일반적 표기)
# 원형(e=0)이면 간단히 v_r = K * sin(2πt/P + φ) 또는 cos-형식. 여기서는 sin 사용.
# K = (2πG/P)^{1/3} * Mp sin i / (Ms+Mp)^{2/3} / sqrt(1-e^2)
K = ((2.0 * math.pi * G / P) ** (1.0 / 3.0)) * (Mp_kg * math.sin(inc)) / ((Ms_kg + Mp_kg) ** (2.0 / 3.0)) / math.sqrt(1.0 - e * e)

# 시간축 설정
T_show = t_window_factor * P
T0 = phase * P  # 시작 위상 이동

t = np.linspace(0.0, T_show, npts)  # 로컬 윈도우

if e == 0.0:
    # 원형 궤도: 단순 사인파
    vr = K * np.sin(2.0 * np.pi * (t + T0) / P)
else:
    # 타원 궤도: 평균이각 M → 편심이각 E → 진근점이각 θ
    M = 2.0 * np.pi * (t + T0) / P  # 평균이각
    # 케플러 방정식: M = E - e sin E → 뉴턴 방법으로 E 풀이
    def solve_E(m):
        E = m  # 초기값
        for _ in range(25):
            f = E - e * np.sin(E) - m
            fp = 1 - e * np.cos(E)
            E -= f / fp
        return E
    E = solve_E(M)
    theta = 2.0 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2.0), np.sqrt(1 - e) * np.cos(E / 2.0))
    vr = K * (np.cos(theta + omega) + e * np.cos(omega))

# ----------------------------- Figures -----------------------------
col1, col2, col3 = st.columns(3, gap="large")

# 1) RV curve
with col1:
    st.subheader("시선속도 곡선 (RV)")
    fig1, ax1 = plt.subplots(figsize=(5, 3))
    ax1.plot(t / DAY, vr)
    ax1.axhline(0, lw=1, ls=":")
    ax1.set_xlabel("시간 (일)")
    ax1.set_ylabel("시선속도 v_r (m/s)")
    ax1.set_title(f"K = {K:,.2f} m/s, P = {P_days:.2f} d")
    fig1.tight_layout()
    st.pyplot(fig1, use_container_width=True)

# 2) Spectrum schematic — 한 줄의 흡수선 이동
with col2:
    st.subheader("흡수선 도플러 편이 (스케치)")
    # 기준 파장 (H-alpha 근처 예시)
    lambda0 = 656.28e-9  # m
    # 중간 시점의 속도(대충 대표)
    v_mid = float(vr[len(vr)//2])
    lambda_shift = lambda0 * (1 + v_mid / C)
    # 간단한 세로선으로 표현
    fig2, ax2 = plt.subplots(figsize=(5, 3))
    # 스펙트럼 바탕 (의미만 전달)
    x = np.linspace(650, 662, 500)  # nm
    y = np.ones_like(x)
    ax2.plot(x, y)
    # 기준선과 이동선
    ax2.axvline(lambda0 * 1e9, ymin=0, ymax=1, ls="--", lw=1)
    ax2.axvline(lambda_shift * 1e9, ymin=0, ymax=1, lw=2)
    ax2.set_xlim(650, 662)
    ax2.set_ylim(0.9, 1.1)
    ax2.set_xlabel("파장 (nm)")
    ax2.set_yticks([])
    direction = "청색편이" if v_mid < 0 else ("적색편이" if v_mid > 0 else "무편이")
    ax2.set_title(f"중간 시점 도플러: {direction} / Δλ ≈ {(lambda_shift - lambda0)*1e12:.3f} pm")
    fig2.tight_layout()
    st.pyplot(fig2, use_container_width=True)

# 3) Orbit view (투영)
with col3:
    st.subheader("궤도 투영 (i에 따른 모양)")
    theta = np.linspace(0, 2*np.pi, 600)
    # 행성 궤도 (항성 기준)
    # 반장축 a_p ≈ a * Ms/(Ms+Mp) 이지만 Mp << Ms 가정으로 시각화는 a만 사용
    x_orb = np.cos(theta)
    y_orb = np.sin(theta) * np.cos(inc)
    fig3, ax3 = plt.subplots(figsize=(5, 3))
    ax3.plot(x_orb, y_orb, lw=2)
    # 현재 위상 위치 마커 (T0를 기준으로)
    th_now = 2.0 * np.pi * ((T0 % P) / P)
    x_now = np.cos(th_now)
    y_now = np.sin(th_now) * np.cos(inc)
    ax3.plot([0], [0], marker="*", markersize=12)  # 항성
    ax3.plot([x_now], [y_now], marker="o")  # 행성
    ax3.set_aspect("equal", adjustable="box")
    ax3.set_title(f"i = {inc_deg}° (원형 투영)")
    ax3.set_xticks([]); ax3.set_yticks([])
    fig3.tight_layout()
    st.pyplot(fig3, use_container_width=True)

# ----------------------------- Notes -----------------------------
st.markdown(
    """
**참고**
- Streamlit Cloud 환경에서는 Tkinter/PyQt 같은 GUI 백엔드를 사용할 수 없습니다. 본 코드는 `matplotlib.use('Agg')`로 설정되어 있으며, 모든 그림은 `st.pyplot`으로 렌더링합니다.
- 원형(e=0) 기본이며, e>0일 때는 표준 RV 공식에 따라 계산합니다. (간단한 수치해로 E를 풉니다.)
- 도플러 스케치의 Δλ는 중간 시점 속도(v_mid)를 사용한 근사입니다.
"""
)