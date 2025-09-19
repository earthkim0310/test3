# Streamlit RV (Radial Velocity) Exoplanet Simulator — Cloud friendly + 재생/일시정지 애니메이션
# - No GUI backends (TkAgg 제거), matplotlib Agg 사용
# - Streamlit 위젯로 제어, st.pyplot으로 렌더링
# - 재생/일시정지, 속도, 리셋 구현 (st.session_state + st_autorefresh)
# - RV 곡선 + 현재 시점 마커, 도플러 선 이동, 궤도 투영과 현재 위치 마커
#
# 요구 패키지 (requirements.txt 예시)
# streamlit
# matplotlib
# numpy
# pandas
# pillow

import math
import time
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
    "Streamlit Cloud 호환 버전 — matplotlib 'Agg' 백엔드 사용, Tkinter/GUI 의존성 제거 + 재생/일시정지 애니메이션"
)

# ----------------------------- Init session_state -----------------------------
if "playing" not in st.session_state:
    st.session_state.playing = False
if "t_start" not in st.session_state:
    st.session_state.t_start = time.time()  # 최근 재생 시작 시간(실시간 now 기준)
if "t_accum" not in st.session_state:
    st.session_state.t_accum = 0.0  # 누적 재생 시간(초)
if "speed" not in st.session_state:
    st.session_state.speed = 1.0

# ----------------------------- Sidebar Controls -----------------------------
with st.sidebar:
    st.header("⚙️ 파라미터")
    Ms = st.number_input("항성 질량 Ms (M☉)", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
    Mp = st.number_input("행성 질량 Mp (M♃)", min_value=0.01, max_value=20.0, value=1.0, step=0.01)
    P_days = st.number_input("공전 주기 P (일)", min_value=0.1, max_value=5000.0, value=365.25, step=0.1)
    inc_deg = st.slider("기울기 i (°)", min_value=0, max_value=90, value=90)
    e = st.slider("이심률 e", min_value=0.0, max_value=0.9, value=0.0, step=0.01)
    omega_deg = st.slider("근일점 인수 ω (°) — e>0일 때 의미", min_value=0, max_value=360, value=90)

    st.markdown("---")
    st.subheader("▶️ 재생 제어")
    speed = st.select_slider("재생 속도", options=[0.25, 0.5, 1.0, 2.0, 4.0], value=float(st.session_state.speed))
    if speed != st.session_state.speed:
        # 속도 변경 시 기준 시각 재정렬
        if st.session_state.playing:
            # 지금까지 경과분 반영 후 속도만 교체
            now = time.time()
            st.session_state.t_accum += (now - st.session_state.t_start) * st.session_state.speed
            st.session_state.t_start = now
        st.session_state.speed = float(speed)

    cols = st.columns(3)
    with cols[0]:
        if st.button("▶︎ 재생" if not st.session_state.playing else "⏸ 일시정지", use_container_width=True):
            if not st.session_state.playing:
                st.session_state.playing = True
                st.session_state.t_start = time.time()
            else:
                # 일시정지: 지금까지 경과분 누적
                now = time.time()
                st.session_state.t_accum += (now - st.session_state.t_start) * st.session_state.speed
                st.session_state.playing = False
    with cols[1]:
        if st.button("⏮ 리셋", use_container_width=True):
            st.session_state.playing = False
            st.session_state.t_accum = 0.0
            st.session_state.t_start = time.time()
    with cols[2]:
        phase0 = st.slider("시작 위상 (0→1)", 0.0, 1.0, 0.0, 0.01)

    st.markdown("---")
    t_window_factor = st.select_slider("시간창 길이 (주기의 배수)", options=[0.5, 1.0, 2.0, 3.0], value=1.0)
    npts = st.slider("샘플 수", min_value=300, max_value=4000, value=1200, step=100)

# ----------------------------- Derived -----------------------------
Ms_kg = Ms * MSUN
Mp_kg = Mp * MJUP
P = P_days * DAY
inc = math.radians(inc_deg)
omega = math.radians(omega_deg)

# 궤도반지름 a (케플러 3법칙)
a = (G * (Ms_kg + Mp_kg) * P * P / (4.0 * math.pi**2)) ** (1.0 / 3.0)  # meters

# RV semi-amplitude K
# 이심률이 1에 가까울 때 수치적 안정성을 위해 제한
e_safe = min(e, 0.999)
K = ((2.0 * math.pi * G / P) ** (1.0 / 3.0)) * (Mp_kg * math.sin(inc)) / ((Ms_kg + Mp_kg) ** (2.0 / 3.0)) / math.sqrt(1.0 - e_safe * e_safe)

# ----------------------------- Time base (animation) -----------------------------
# 누적 경과시간(초): 일시정지/재생/속도를 반영
if st.session_state.playing:
    now = time.time()
    t_elapsed = st.session_state.t_accum + (now - st.session_state.t_start) * st.session_state.speed
else:
    t_elapsed = st.session_state.t_accum

# 시작 위상(phase0)을 절대시간으로 변환
# phase0는 0~1 범위의 위상이므로 주기에 곱해서 시간으로 변환
T0 = phase0 * P

# 현재 절대시간 (시뮬레이터 내부 기준)
T_curr = T0 + t_elapsed

# 위상이 주기 범위를 벗어나지 않도록 정규화
T_curr = T_curr % P

# 시각화용 시간창
T_show = t_window_factor * P
# 현재 시점을 창 중앙에 두고 좌우로 펼침
left = T_curr - 0.5 * T_show
right = T_curr + 0.5 * T_show

# 샘플링 시간축
t = np.linspace(left, right, npts)

# ----------------------------- Kepler solver -----------------------------
def solve_E(m, e):
    # 벡터 m 입력을 허용하는 뉴턴 방법 구현
    m = np.asarray(m)
    E = m.copy()
    
    # 수치적 안정성을 위한 개선된 초기값
    if e < 0.8:
        E = m + e * np.sin(m) / (1.0 - e * np.cos(m))
    else:
        E = m + e * np.sin(m)
    
    # 뉴턴 방법으로 수렴
    for _ in range(50):  # 반복 횟수 증가
        f = E - e * np.sin(E) - m
        fp = 1.0 - e * np.cos(E)
        
        # 0으로 나누기 방지
        delta = np.where(np.abs(fp) > 1e-12, f / fp, 0.0)
        E -= delta
        
        # 수렴 조건 확인
        if np.all(np.abs(delta) < 1e-10):
            break
    
    return E

# 평균이각
M = 2.0 * np.pi * (t / P)

if e == 0.0:
    # 원형 궤도: 사인파로 단순화
    vr = K * np.sin(2.0 * np.pi * (t / P))
    # 현재 시점 속도
    v_now = float(K * np.sin(2.0 * np.pi * (T_curr / P)))
    # 현재 위상(궤도 표시용)
    theta_now = 2.0 * np.pi * (T_curr / P)  # [0, 2π] 주기적 해석
else:
    # 타원 궤도: 표준 RV 식 사용
    try:
        E = solve_E(M, e_safe)
        theta = 2.0 * np.arctan2(np.sqrt(1 + e_safe) * np.sin(E / 2.0), np.sqrt(1 - e_safe) * np.cos(E / 2.0))
        vr = K * (np.cos(theta + omega) + e_safe * np.cos(omega))

        # 현재 시점
        M_now = 2.0 * np.pi * (T_curr / P)
        E_now = solve_E(M_now, e_safe)
        theta_now = 2.0 * np.arctan2(np.sqrt(1 + e_safe) * np.sin(E_now / 2.0), np.sqrt(1 - e_safe) * np.cos(E_now / 2.0))
        v_now = float(K * (np.cos(theta_now + omega) + e_safe * np.cos(omega)))
    except (ValueError, ZeroDivisionError, OverflowError):
        # 수치 오류 발생 시 원형 궤도로 대체
        vr = K * np.sin(2.0 * np.pi * (t / P))
        v_now = float(K * np.sin(2.0 * np.pi * (T_curr / P)))
        theta_now = 2.0 * np.pi * (T_curr / P)

# ----------------------------- Figures -----------------------------
# 애니메이션을 위한 빈 컨테이너 생성
plot_container = st.empty()

with plot_container.container():
    col1, col2, col3 = st.columns(3, gap="large")

    # 1) RV curve + current-time marker
    with col1:
        st.subheader("시선속도 곡선 (RV)")
        fig1, ax1 = plt.subplots(figsize=(5.5, 3.3))
        ax1.plot((t - T_curr) / DAY, vr)  # 현재 시점을 0으로 이동한 로컬 좌표(일)
        ax1.axvline(0, lw=1.2, ls=":")
        ax1.axhline(0, lw=1.0, ls=":")
        ax1.set_xlabel("현재시점 기준 시간 (일)")
        ax1.set_ylabel("시선속도 v_r (m/s)")
        ax1.set_title(f"K = {K:,.2f} m/s, P = {P_days:.2f} d")
        fig1.tight_layout()
        st.pyplot(fig1, use_container_width=True)

    # 2) Spectrum schematic — 한 줄의 흡수선 이동 (현재 시점 속도 기반)
    with col2:
        st.subheader("흡수선 도플러 편이 (현재 시점)")
        lambda0 = 656.28e-9  # m (H-alpha 근처)
        lambda_shift = lambda0 * (1.0 + v_now / C)
        fig2, ax2 = plt.subplots(figsize=(5.5, 3.3))
        x = np.linspace(650, 662, 500)  # nm
        y = np.ones_like(x)
        ax2.plot(x, y)
        ax2.axvline(lambda0 * 1e9, ymin=0, ymax=1, ls="--", lw=1)  # 기준선
        ax2.axvline(lambda_shift * 1e9, ymin=0, ymax=1, lw=2)       # 이동선
        ax2.set_xlim(650, 662)
        ax2.set_ylim(0.9, 1.1)
        ax2.set_xlabel("파장 (nm)")
        ax2.set_yticks([])
        direction = "청색편이" if v_now < 0 else ("적색편이" if v_now > 0 else "무편이")
        ax2.set_title(f"v_now = {v_now:,.2f} m/s → {direction}\nΔλ ≈ {(lambda_shift - lambda0)*1e12:.3f} pm")
        fig2.tight_layout()
        st.pyplot(fig2, use_container_width=True)

    # 3) Orbit view (투영) + 현재 위치 마커
    with col3:
        st.subheader("궤도 투영 (i에 따른 모양)")
        th = np.linspace(0, 2*np.pi, 600)
        x_orb = np.cos(th)
        y_orb = np.sin(th) * np.cos(inc)
        fig3, ax3 = plt.subplots(figsize=(5.5, 3.3))
        ax3.plot(x_orb, y_orb, lw=2)
        x_now = np.cos(theta_now)
        y_now = np.sin(theta_now) * np.cos(inc)
        ax3.plot([0], [0], marker="*", markersize=12)  # 항성
        ax3.plot([x_now], [y_now], marker="o")  # 행성 현재 위치
        ax3.set_aspect("equal", adjustable="box")
        ax3.set_title(f"i = {inc_deg}° (원형/타원 투영)")
        ax3.set_xticks([]); ax3.set_yticks([])
        fig3.tight_layout()
        st.pyplot(fig3, use_container_width=True)

# ----------------------------- Autorefresh when playing -----------------------------
# 재생 중일 때만 주기적 리렌더링 (100ms)
if st.session_state.playing:
    time.sleep(0.1)  # 100ms 대기
    # plot_container를 사용하여 애니메이션 효과 구현
    plot_container.empty()
    # ScriptRunContext 경고는 무시해도 됨 (기능에 영향 없음)
    st.rerun()

# ----------------------------- Notes -----------------------------
st.markdown(
    """
**사용 팁**
- ▶︎ 재생/⏸ 일시정지, ⏮ 리셋 버튼으로 시간 흐름을 제어함. 재생 속도를 0.25×~4×로 변경 가능함.
- 좌측 플롯은 현재 시점을 가운데(0일)로 한 로컬 시간 창을 보여주며, 세로 점선이 현재 시점을 나타냄.
- 가운데 플롯은 현재 시점 속도(v_now)에 따른 흡수선의 도플러 이동을 시각화함.
- 우측 플롯은 i(기울기)에 따른 궤도 투영과 현재 행성 위치를 표시함.
- 이심률 e>0일 때는 표준 RV 공식 v_r = K[cos(θ+ω) + e cos ω]을 사용하며, 케플러 방정식을 수치해로 풀이함.
"""
)