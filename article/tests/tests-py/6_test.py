"""
Om(z)-диагностика с модельно-независимой нормировкой H₀.

В предыдущей версии H₀ извлекался подгонкой хронометров формой ΛCDM,
что вносит зависимость от модели в нормировку. Здесь H₀ задаётся тремя
независимыми способами, и проверяется устойчивость вывода:

  (1) непараметрическая экстраполяция H(z)→0 (линейная, без модели);
  (2) внешний H₀ = 73.0 (SH0ES, локальная лестница расстояний);
  (3) внешний H₀ = 67.4 (Planck, ранняя Вселенная).

Если знак Δχ²(DSIC−ΛCDM) одинаков при всех нормировках, вывод устойчив.
Если знак меняется, Om(z)-тест не разрешает модели на текущих данных.

Зависимости: numpy, scipy, matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

c = 299792.458

# ---- Данные: космические хронометры H(z) (компиляция Moresco) ----
CC = np.array([
 [0.07,69.0,19.6],[0.09,69.0,12.0],[0.12,68.6,26.2],[0.17,83.0,8.0],
 [0.179,75.0,4.0],[0.199,75.0,5.0],[0.20,72.9,29.6],[0.27,77.0,14.0],
 [0.28,88.8,36.6],[0.352,83.0,14.0],[0.38,83.0,13.5],[0.40,95.0,17.0],
 [0.4004,77.0,10.2],[0.425,87.1,11.2],[0.445,92.8,12.9],[0.47,89.0,49.6],
 [0.4783,80.9,9.0],[0.48,97.0,62.0],[0.593,104.0,13.0],[0.68,92.0,8.0],
 [0.75,98.8,33.6],[0.781,105.0,12.0],[0.80,113.1,28.5],[0.875,125.0,17.0],
 [0.88,90.0,40.0],[0.90,117.0,23.0],[1.037,154.0,20.0],[1.30,168.0,17.0],
 [1.363,160.0,33.6],[1.43,177.0,18.0],[1.53,140.0,14.0],[1.75,202.0,40.0],
 [1.965,186.5,50.4],
])
z_cc,H_cc,sH_cc = CC[:,0],CC[:,1],CC[:,2]

# ---- BAO DESI DR2: D_H/r_d ----
BAO=np.array([[0.510,21.863,0.425],[0.706,19.455,0.330],[0.930,17.641,0.193],
              [1.317,14.176,0.221],[1.491,12.817,0.516],[2.330,8.632,0.101]])
z_b,DHrd,sDHrd=BAO[:,0],BAO[:,1],BAO[:,2]
rd=147.1

# ---- Модели ----
mu0_fit,Om_fit=0.7615,0.3044
def b(m): return np.sqrt(1-m*m)/(m*m)
def mu_e(zv,mu0): return brentq(lambda m:b(m)-b(mu0)/(1+zv),1e-12,1-1e-15)
def E_DSIC(zv):
    me=mu_e(zv,mu0_fit); return ((2-me**2)/(1-me**2))/((2-mu0_fit**2)/(1-mu0_fit**2))
def E_LCDM(zv): return np.sqrt(Om_fit*(1+zv)**3+(1-Om_fit))
def Om_of_E(zv,E): return (E**2-1)/((1+zv)**3-1)

# ---- (1) Непараметрическая экстраполяция H0: H(z)=H0+a·z на z<0.5 ----
m=z_cc<0.5
W=1/sH_cc[m]**2; A=np.vstack([np.ones(m.sum()),z_cc[m]]).T
H0_extrap=np.linalg.solve(A.T@(A*W[:,None]), A.T@(W*H_cc[m]))[0]
print(f"(1) H₀, непараметрическая экстраполяция = {H0_extrap:.1f}")
print(f"(2) H₀, внешний SH0ES = 73.0")
print(f"(3) H₀, внешний Planck = 67.4")

# ---- При заданном H0 строит Om(z) из данных и вычисляет Δχ² ----
def analyze(H0, label):
    # Om(z) из хронометров
    E_d=H_cc/H0; Om_d=Om_of_E(z_cc,E_d)
    dOm=np.abs(2*E_d/((1+z_cc)**3-1))*(sH_cc/H0)
    # Om(z) из BAO
    H_bao=c/(DHrd*rd); E_b=H_bao/H0; Om_b=Om_of_E(z_b,E_b)
    dOm_b=np.abs(2*E_b/((1+z_b)**3-1))*( (H_bao*(sDHrd/DHrd))/H0 )
    # χ² каждой кривой
    def chi2(Ef):
        c1=np.sum(((Om_d-np.array([Om_of_E(z,Ef(z)) for z in z_cc]))/dOm)**2)
        c2=np.sum(((Om_b-np.array([Om_of_E(z,Ef(z)) for z in z_b]))/dOm_b)**2)
        return c1+c2
    cD,cL=chi2(E_DSIC),chi2(E_LCDM)
    verdict=('ближе DSIC' if cL-cD>2 else 'ближе ΛCDM' if cD-cL>2 else 'не различить')
    print(f"  {label} (H₀={H0:.1f}):  χ²_DSIC={cD:.1f}  χ²_ΛCDM={cL:.1f}  Δχ²={cD-cL:+.1f}  → {verdict}")
    return (z_cc,Om_d,dOm,z_b,Om_b,dOm_b,cD-cL)

print("\nУСТОЙЧИВОСТЬ ВЫВОДА К ВЫБОРУ H₀:")
res_extrap=analyze(H0_extrap,"экстраполяция")
res_sh0es=analyze(73.0,"SH0ES     ")
res_planck=analyze(67.4,"Planck    ")

# ---- Скан по H0: знак Δχ² как функция H0 ----
print("\nСКАН Δχ²(DSIC−ΛCDM) ПО H₀ (устойчивость знака):")
H0scan=np.linspace(64,76,13)
signs=[]
for H0 in H0scan:
    E_d=H_cc/H0; Om_d=Om_of_E(z_cc,E_d); dOm=np.abs(2*E_d/((1+z_cc)**3-1))*(sH_cc/H0)
    H_bao=c/(DHrd*rd); E_b=H_bao/H0; Om_b=Om_of_E(z_b,E_b)
    dOm_b=np.abs(2*E_b/((1+z_b)**3-1))*((H_bao*(sDHrd/DHrd))/H0)
    def chi2(Ef):
        return (np.sum(((Om_d-np.array([Om_of_E(z,Ef(z)) for z in z_cc]))/dOm)**2)+
                np.sum(((Om_b-np.array([Om_of_E(z,Ef(z)) for z in z_b]))/dOm_b)**2))
    dchi=chi2(E_DSIC)-chi2(E_LCDM); signs.append(dchi)
    print(f"  H₀={H0:.0f}: Δχ²={dchi:+.1f}")
signs=np.array(signs)
allneg=np.all(signs<0); allpos=np.all(signs>0)
print("\nВЕРДИКТ:")
if allneg:
    print("  Δχ²<0 при всех H₀ → данные устойчиво ближе к DSIC")
    print("  (слабо, но без артефакта нормировки)")
elif allpos:
    print("  Δχ²>0 при всех H₀ → данные устойчиво ближе к ΛCDM")
else:
    print("  знак Δχ² меняется с H₀ → Om(z)-тест не разрешает модели")
    print("  на текущих данных")
print(f"  (диапазон Δχ²: {signs.min():+.1f} … {signs.max():+.1f})")

# ---- График с тремя нормировками ----
zgrid=np.linspace(0.05,2.5,200)
plt.figure(figsize=(10,6))
plt.plot(zgrid,[Om_of_E(z,E_LCDM(z)) for z in zgrid],color='steelblue',lw=2,label='ΛCDM (плоская)')
plt.plot(zgrid,[Om_of_E(z,E_DSIC(z)) for z in zgrid],color='crimson',lw=2,label='DSIC (U-образная)')
for (H0,col,mk) in [(H0_extrap,'black','o'),(73.0,'gray','^'),(67.4,'darkgreen','v')]:
    E_d=H_cc/H0; Om_d=Om_of_E(z_cc,E_d)
    plt.scatter(z_cc,Om_d,c=col,s=22,alpha=0.5,marker=mk,label=f'хронометры (H₀={H0:.0f})')
plt.xlabel('z'); plt.ylabel('Om(z)'); plt.ylim(-0.1,0.9); plt.grid(alpha=0.3)
plt.title('Om(z): устойчивость к нормировке H₀'); plt.legend(fontsize=8)
plt.tight_layout(); plt.savefig('Om_z_robust.png',dpi=130); plt.show()
print("\nГрафик сохранён: Om_z_robust.png")
