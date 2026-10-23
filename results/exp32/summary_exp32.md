# Exp32 — Tutmayan modellerin ölçümle düzeltilmesi ve SKY130'un kendi 4 K ölçümü

Tarih: 2026-10-16. Kod: `experiments/exp32_corrections/`
- `corr32.py`: C1, C2, C3, C5;
- `sub32.py`: alt bant eşiğinin katkıya göre değişimi;
- `nist32.py`: SKY130 4 K verisi ve yeniden kalibrasyon.

Öngörüler hesaptan önce yazıldı:
- `predictions_exp32.json` (15:26);
- `predictions_exp32_nist.json` (C4'ün işlemsel tanımları, analizden önce).

Sonuçlar:
- `corrections.json`;
- `nist_sky130.json`;
- `subband_NA*.json`;
- `nist_sky130/` (kullanılan ham ölçüm dosyaları: 4 K ve 300 K).

Şekil: `fig_exp32_nist.png`.

**Neden.** Stop kancası şunlara itiraz etti:
- SKY130'a özgü hücre sabitleri ölçülmemişti;
- hücre terimi başka bir aygıt sınıfıyla sınanmıştı;
- 7 sonuç karışık, 5'i tutmuyordu.

Kural: ölçümle çelişen model düzeltilir, düzeltilmiş model düzeltmede kullanılmayan bir ölçümle yeniden sınanır.

## Kısa sonuç

- **SKY130'un kendi 4 K ölçümü bulundu ve kullanıldı.** NIST mds2-2997 (Akturk, Hoskins, Shrestha vd.): ETest yongası, Lakeshore kapalı
  devre kriyojenik prob istasyonu, örnek ve çevresi **4 K**, Agilent 4156B. Oda sıcaklığı karşılığı aynı aygıt kimliğiyle (8392_7_8)
  google/skywater-pdk-sky130-raw-data'dan.
  - **C4a ✓:** kink kriyojenik. nfet_01v8 W7/L0,15'te kink oranı 4 K'de 2,77, 300 K'de 1,15. W7/L8'de 32,7 ve 1,25. 4 K'deki 10
    aygıtın hepsinde kink var. Mekanizma Exp12'nin mekanizması: I_D artışı I_B'nin ~2000 katı, yani delikler donmuş gövdeyi yükseltip
    eşiği düşürüyor.
  - **C4b ✗ (kayıtlı ölçüte göre):** eşik üstü kink 1,25–1,525 V'ta başlıyor, 0,6–1,2 V aralığında değil. Gövde akımı 1,375–1,8 V'ta
    başlıyor. Not: bu ölçüt eşik üstü kinki, eşik altı kilitlenmenin kenet gerilimiyle karıştırıyordu; ikisi ayrı büyüklükler.
  - **C4c — SKY130 sabitleri ölçülmüş girdilerle yeniden hesaplandı.**
    - Exp17'nin çoğaltma varsayımı (300 K BSIM4 × 10) ölçümün 6,6 (0,8 V) ile 22 (1,8 V) katıydı.
    - Ölçülmüş M ve ölçülmüş 4 K transfer eğrisiyle: **V_M = 0,701 V** (yayımlanan 0,802), **j_on = 64 µA/µm** (yayımlanan 2,57).
    - Hücreler daha dar (w_s 20–54 µm) ve her biri 1,0–2,6 mA taşıyor.
    - Benim öngördüğüm "V_M 1,2–1,6 V" tutmadı.
  - **Tutarlılık ✓:** düzeltilmiş modelde düşük durum ≥ 2,3 V'a kadar kararlı (V_up). Ölçümde eşiğin 0,2 V altında akım 1,8 V'a
    kadar pA düzeyinde kalıyor; kendiliğinden kilitlenme yok. Doğum için tohum gerekir (Exp18).
  - **Histerezis:** ölçüm yok (yalnız ileri tarama).
- **C1 ✓ (Exp0–2 düzeltildi).** Zayıf kapılı adada doğru yasa öz-kapasitans C = 2π ε₀ ε_ox d.
  - Shin'in TEM boyutuyla (2 nm) 0,434 aF öngörüyor; ölçülen 0,42 aF (+%3,3).
  - 5 nm adada e²/C 300 K'de 5,7 kT (ölçüm: tek tepe), 200 K'de 8,6 kT (ölçüm: belirgin blokaj).
  - Ölçülen kapı kaldıracı 0,22, Exp1 geometrisinin 0,76'sı.
  - Düzeltilmiş eşikler: 300 K ≤ 1,4 nm, 77 K ≤ 5,6 nm, 4 K ≤ 107 nm.
- **C2 ✓ (Exp6 düzeltildi).** Ölçülmüş Pelgrom kuralı (A_VT 4,15–5 mV·µm) ile bağımsız SKY130 PDK katsayısı (3,36 mV·µm) 1,24–1,49 kat
  içinde. +30 mV tohum yine 7,3–10,9σ.
- **C3 ✗ (Exp22 kalıyor).** Kendi Poisson–Schrödinger çözücümüz 3·10¹⁵ – 1,7·10¹⁷ katkının hepsinde tek alt bandı ≥ 1,4·10¹³ cm⁻²'ye
  kadar sürdürüyor; ölçülen eşik ~6·10¹². Neden katkı farkı değil: çözücü (Hartree, sert duvar, değişim-korelasyon yok) vadiler arası
  aralığı ~2 kat büyük veriyor. Ölçülen değer kullanılır; itki yolunda yok.
- **C5 (CL2 yönü).** SKY130 modelinde j_high, V_D ile artıyor; model "hücre sayısı artınca j artar" diyor. ESD ölçümünde (ısıl, 300 K)
  azalıyor. Aygıt sınıfı farkı; SKY130 için kayıtlı öngörü.

## SKY130 4 K ölçüm özeti (nfet_01v8, V_B = 0)

| aygıt | V_T (V_D 0,1 V) | en küçük alt eşik eğimi | kink oranı (V_G 0,9 V) | kink başlangıcı | I_B > 1 nA |
|---|---|---|---|---|---|
| W7/L0,15 | 0,797 V | 5,8 mV/dek | 2,77 | 1,275 V | 1,425 V |
| W1,68/L0,15 | 0,838 V | 11,6 mV/dek | 3,48 | 1,25 V | 1,525 V |
| W0,42/L0,15 | 0,840 V | 7,4 mV/dek | 3,99 | 1,25 V | 1,65 V |
| W7/L8 | 0,687 V | 9,7 mV/dek | 32,7 | 1,425 V | 1,8 V |
| lvt W7/L0,15 | 0,709 V | 4,8 mV/dek | 1,53 | 1,375 V | 1,375 V |

300 K'de aynı W7/L0,15 aygıtında kink oranı 1,15 (V_G 0,72 V) ve 0,73 (1,08 V); W7/L8'de 1,25 ve 0,30. Eğim değerleri 0,01 V'luk tarama
adımıyla sınırlı.

## Ölçülmüş girdilerle hücre sabitleri (Exp17 yasası)

| | V_M | j_on | λ_c (r_D 5–20) | w_s (λ 40–120 µm) | I_c | hücre başına itki (In⁺, 5 kV, γ 0,9) |
|---|---|---|---|---|---|---|
| Exp17 (yayımlanan; varsayımlı 4 K düğmeleri) | 0,802 V | 2,57 µA/µm | 25–50 µm | 56–219 µm | 0,11–0,42 mA | 11–41 µN |
| SKY130 4 K ölçülmüş M + transfer eğrisi | 0,701 V | 64 µA/µm | 9,8–19,7 µm | 20,5–54 µm | 0,99–2,6 mA | 97–255 µN |

Hâlâ varsayım olan: ileri gövde ön geriliminde gövde etkisi (BSIM4) ve gövde diyot sızıntısı (V_on 1,05 V, V_s 5 mV). Kilitlenmenin
(filamanın) kendisi geniş bir SKY130 aygıtında doğrudan ölçülmüş değil; NIST verisinde en geniş nfet_01v8 7 µm ve balast yok.

## Bütün verinin son durumu (Exp30 + Exp31 + Exp32)

| sonuç | sayı | deneyler |
|---|---|---|
| ✓ (ölçümle düzeltilenler dahil) | 16 | Exp0, 1, 2 (C1); 6 (C2); 3, 4, 5, 8, 11, 13, 15, 18, 19, 21, 23, 27 |
| karışık / kısmen | 6 | Exp10, 12 (C4a ✓ / C4b ✗), 14, 17 (3 ✓ / CL2 ✗; sabitler ölçülmüş girdilerle), 20, 24 |
| ✗ | 3 | Exp22 (alt bant çözücüsü), 25 (GHZ), 26 (eşevrelilik) |
| yöntem | 3 | Exp7, 9, 16 |

Üç ✗'in hiçbiri itki denklemine girmiyor; kullanılırlarsa ölçülen değerleri kullanılır.

## Doğru yol ve denklem (bütün terimler ölçülmüş girdilerle)

  F = γ · N · I_c · √(2 m V / q),   I_c ≈ 0,75 · j_on · w_s

| terim | dayanak |
|---|---|
| γ, √(2mV/q) | NEXT 40 nokta (%1,45); In FEEP (%1,6–4,3) |
| j_on = 64 µA/µm, V_M = 0,70 V | SKY130 4 K ölçülmüş M ve transfer eğrisi |
| w_s, N(I) | Exp17 yasası; yapısı silisyum filaman ölçümleriyle doğrulandı (Exp31) |
| yayıcı | alanla açılan atom yayımı (FEEP); programlanabilir SiOx yayıcı |
| kriyojenik çalışma | SKY130'da kink 4 K'de var, 300 K'de yok; 180 nm'de 30 K'de var, 40 K'de yok |

## Kaynaklar

- NIST mds2-2997, *CMOS Sky130 Primitives measured at cryogenic temperatures*, sürüm 2.0.0 (Akturk, Li, Saligane, Riem, Adam,
  Hoskins, Shrestha), https://doi.org/10.18434/mds2-2997. Dosyalar https://data.nist.gov/cch/... önbellek adresinden indirildi; ETest
  4 K, 23–24 Mart 2022.
- google/skywater-pdk-sky130-raw-data (oda sıcaklığı, aynı aygıt kimliği),
  https://github.com/google/skywater-pdk-sky130-raw-data
- Shin vd., arXiv 1001.3724 (C1)
- Pelgrom kuralı (Berkeley EE240B); SKY130 PDK vth0_slope (C2)
- Ray ve Martin, IBM (C3)
- Ilani vd., PRL 84, 3133 (2000) ve Dultz–Jiang, PRL 84, 4689 (2000): sıkıştırılabilirlik ölçümleri GaAs'ta; Si için eşik sayısı yok
  (Exp10 kısmen)
