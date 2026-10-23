# Exp22 — Elimizdeki veriyle silikonda (SKY130) kuantum durumu

Tarih: 2026-10-15. Kod: `experiments/exp22_quantum/q22.py` (alt bantlar, Landau, yüklenme, hücre sayımı),
`experiments/exp22_quantum/cv22.py` (ngspice C–V, klasik/kuantum C–V, tünelleme, fig). Şekil: `fig_exp22.png`.

## Dürüst çerçeve

Elimizde ölçüm verisi yok. Var olanlar: resmî SKY130 model kartı (SkyWater'ın silikon ölçümlerine uydurulmuş), projenin
çözücüleri (Poisson–Schrödinger, FEM) ve önceki deneyler. Bu yüzden "gözlem" burada iki anlama geliyor:

- dökümhane verisinde (kartta) bir kuantum izi var mı;
- hangi ölçümle doğrudan görülür.

Kübit düzeyinde (eşevreli, kontrol edilen) bir kuantum durumu bu veriyle gösterilemez. Önceki deneylerin (Exp0–10)
sonucu değişmedi: oda sıcaklığında SKY130'da kübit yok.

## Bulgular

1. **İnversiyon katmanı bir kuantum durumudur (her SKY130 transistöründe, oda sıcaklığında da).**
   - Öz-uyumlu Poisson–Schrödinger (t_ox 4,148 nm, N_A 1,7e17, 300 K) sonuçları:
     - E₁ − E₀ = 43–91 meV (n_s 6e10–9,5e12 cm⁻²);
     - 300 K'de elektronların %47–72'si en düşük alt bantta;
     - 4 K'de ölçülen bütün aralıkta %100 tek alt bant (elektrik kuantum limiti, E_F − E₀ ≤ 60 meV < E₁ − E₀).
   - Yük merkezi arayüzden 1,2–1,9 nm uzakta (klasik katman arayüze yapışır).
   - Sınırlar: sonsuz arayüz bariyeri, çok-cisim (değiş-tokuş) etkisi yok. Bu ikisi merkezi biraz abartır.
2. **Dökümhane kartında iz (model düzeyinde):** SKY130 kartı capMod = 2 kullanıyor. Bu, BSIM4'ün kuantum
   "yük kalınlığı" (karanlık bölge) modeli; TOXE = TOXP = 4,148 nm. ngspice ile kartın C–V'si (W = 10 µm,
   L = 50 − 25 µm farkı, 100 kHz), bükülme noktasının 0,9 V üstünde (C_ox = 0,8325 µF/cm²):

   | eğri | C/C_ox |
   |---|---|
   | resmî kart | 0,926 |
   | klasik (kuantumsuz) | 0,950 |
   | kuantum P-S | 0,875 |

   - Kartın inversiyon kapasitesi klasiğin %2,4 altında, yani kuantum hapsi yönünde.
   - Ama bu bir ölçüm değil. BSIM4'ün varsayılan katsayılı (ADOS = BDOS = 1) kuantum düzeltmesi; dökümhanenin ölçümle
     ne kadar uydurduğunu bilmiyoruz.
   - Kapı tünellemesi ve GIDL kartta kapalı (igcmod = igbmod = 0, agidl = 0). Dökümhane bu kuantum etkilerini modellememiş.
3. **Tünelleme, oda sıcaklığında ölçülebilir bir kuantum etkisi:**
   - 4,148 nm SiO₂'den doğrudan tünelleme (Schuegraf–Hu WKB, bariyer 3,1 eV) 1,8 V oksit geriliminde J ≈ 5e-12 – 5e-10
     A/cm² (m_ox 0,5–0,4 m₀; belirsizlik ±1–2 mertebe).
   - Isıl emisyon ~1e-45 A/cm². Oksitten geçen akımın tamamı tünelleme (ya da tuzak destekli tünelleme).
   - 1 mm²'lik kapasitörde 46 fA – 4,9 pA: elektrometreyle ölçülebilir.
   - Parmak izi: 300 K ↔ 77 K ↔ 4 K arasında akım neredeyse değişmez; ısıl yollar ise mertebelerce düşer.
4. **4 K + manyetik alan (2B elektron gazı):**
   - ħω_c = 0,61 meV/T; kT(4 K) = 0,34 meV, yani B > 0,57 T yeter.
   - Ama Shubnikov–de Haas salınımı için μB > 1 gerekir: μ = 1000 cm²/Vs'de 10 T, 2000'de 5 T.
   - Dolum faktörü ν = 4 alanı: n_s = 1e12'de 10,3 T.
   - Ağır katkılı kanalda 4 K hareketliliği bilinmiyor; 9–14 T'lik bir laboratuvar mıknatısı gerekir.
5. **Coulomb blokajı (önceki FEM sabiti, E_C/k_B = 585–1171 K·nm / a):**
   - SKY130'un 150 nm'lik kapı-tanımlı adası: E_C/k_B = 3,9–7,8 K. Net blokaj (E_C > 10 kT) için ≤ 0,4–0,8 K, yani
     seyreltme soğutucusu gerekir.
   - 10 nm'lik ada 59–117 K verirdi, ama standart süreçte yapılamaz.
6. **Akım filamanı hücreleri klasik nesneler:**
   - 100 µm'lik bir hücrede kanalda ~1,5e5 elektron, gövdede ~3e5 delik.
   - Hücre süperpozisyonu ya da hücreler arası eşevreli bağlaşım yok.

## Sonuç

- SKY130 silikonunda kuantum durumları var: inversiyon alt bantları, tünelleme, 4 K'de tek alt bantlı 2B gaz ve
  ≤ 1 K'de tek elektron adaları. Elimizdeki veride ise yalnızca model düzeyinde bir iz var (capMod 2 karanlık bölge,
  klasiğin %2,4 altı).
- Doğrudan gözlem için en ucuz ölçüm (oda sıcaklığı): büyük alanlı SKY130 kapasitöründe yüksek çözünürlüklü C–V
  (%1'lik hassasiyet yeter) ve elektrometreyle kapı tünelleme akımının sıcaklıktan bağımsızlığı.
- 4 K: elektrik kuantum limiti (C–V'de ikinci alt bant basamağı yok, n_s ≤ 9,5e12) ve ≥ 9 T'de SdH.
- ≤ 1 K: Coulomb salınımı (Aşama 4).
- Bunlar öngörü; hiçbiri burada ölçülmedi.
