# Exp20 — Hücre sayısını darbe/çukurla kontrol, halka cihaz, ortak düğüm kapasitesi ve çözücü düzeltmesi

Tarih: 2026-10-15. Kod: `experiments/exp20_pulse_control/run20.py`, `recheck17.py`; çözücü: `experiments/exp12_filament/spatial.py`
(`periodic`, `C_c`, `implicit_load`). Tahmin dosyaları koşulardan önce yazıldı: `predictions_*.json`.

## Kısa sonuç

- Halka (kenarsız) cihazda yavaş taramada en güçlü tohum kazanıyor (19/20). Exp19'daki "dış tohum" yanlılığı geometridendi.
- Ani basamakta kaç hücre doğacağını basamak yüksekliği belirliyor (1,05–1,1 V: 2; 1,2 V: 3; 1,4 V: 4). Kazananlar en güçlü
  k tohum değil (tek hücre en büyük sürüşüyle 12/20, 9/20, 5/19). Sonradan bulunan "eşik altı toplam ağırlık" öngörücüsü
  daha iyi, ama yeni alanlarda (ring2) önceden kaydedilen ≥ %70 eşiğini 3 basamaktan 2'sinde tuttu.
- İki düğme: darbe yüksekliği doğumları, tutma gerilimi hayatta kalanları ("taşıma kapasitesi") belirliyor.
- Çip üstü ortak düğüm kapasiteleri (R_L·C ≤ 5 ms) sonucu değiştirmiyor; anlık yük doğrusu varsayımı geçerli.
- Açık (explicit) çözücüde bir yapaylık bulundu ve düzeltildi (örtük yük doğrusu). Exp17'nin yayımlanan sayıları yeni
  çözücüyle büyük ölçüde korunuyor (bkz. "Exp17 yeniden sınaması"); bir koşu ve bir geçerlilik-dışı cümle değişiyor.
- Ani zorla küçültme (kapı çukurunun bırakılması ya da ≥ 0,2 V V_DD düşüşü) hücreleri eşzamanlı bölüyor
  (5,0 V, −60 mV: 3 → 7).

## A — Halka cihaz (alanlar 500–519; tohumlar 75/225/375/525 µm, 45–55 mV)

| koşul | k (hücre) | tohum dışı hücre | kazanan = tek hücre en büyük sürüşle ilk k | kazanan = toplam ağırlıkla ilk k |
|---|---|---|---|---|
| yavaş tarama | 1 (20/20) | 0 | 19/20 (A1 ✓) | 15/20 |
| basamak 1,05 V | 2 (20/20) | 0 | 12/20 | 15/20 |
| basamak 1,1 V | 2 (20/20) | 0 | 9/20 | 14/20 |
| basamak 1,2 V | 3 (19/20) | 0 | 5/19 | 14/19 |
| basamak 1,4 V | 4 (20/20) | 0 | — | — |

- A1 ✓ (yavaş taramada en güçlü tohum). A2 ✗ (basamakta kazananlar en güçlü k tohum değil).
- Toplam ağırlık Σexp(sürüş/nφ) sonradan bulundu; yeni alanlarda (520–539, `predictions_ring2.json`, eşik ≥ %70)
  15/20 ✓, 13/20 = %65 ✗, 15/19 ✓. Her basamakta tek hücre öngörücüsünden daha iyi (1,2 V: 15/19'a karşı 7/19).

## B — Ortak düğüm kapasitesi C_c

- R_L·C ≤ 5 ms: k ve ani öldürme doğumları C_c = 0 ile aynı (B1 ✓, B2 ✓).
- R_L·C ≥ 0,5 s: basamak yumuşuyor, k ortalaması 1,1 V'ta 1,2 ve 1,4 V'ta 1,8, ani öldürme doğumu 0.
- Sonuç: çip üstü kapasiteler için anlık yük doğrusu geçerli.

## C — Darbe ve "taşıma" (carry)

- Darbe deneyinde tasarım hatası: V_p = 1,1 V için tutma gerilimi (1,2 V) darbeden yüksekti, protokol merdivene dönüştü.
  Öğrenilen: yarış yaklaşık 50 ms içinde belirleniyor.
- Carry (0,1 s 1,4 V, sonra V_hold): hayatta kalan hücre sayısı tutma gerilimiyle artıyor (C3 ✓).

| V_hold (V) | 0,95 | 1,0 | 1,1 | 1,2 | 1,3 | 1,4 |
|---|---|---|---|---|---|---|
| hayatta kalan | 1 | 1 | 2 | 3,6 | 4 | 4 |

## D — Çukurlar (dip) ve çözücü yapaylığı

**Açık çözücüyle D koşuları geçersiz.** 54 koşunun hepsi 3 → 3 çıktı. Teşhis: ani akım fazlasından sonra
(kapı çukurunun bırakılması) donmuş-V_ds yinelemesi denge olmayan bir durumda "yakınsıyordu" (V_b ≈ 0,7, iç V_D 0,74 V; gerçek
ray 1,23 V, artık +1,5e-10). Tarama: yapaylık yalnızca çukur koşularında (12/54). Exp19 basamak, halka, kapasite, carry, darbe
ve Exp18 yenilenme koşuları temiz.

**Düzeltme:** `implicit_load=True`: ortak düğüm Newton adımına örtük girer (rank-1 Sherman–Morrison, 2B çözücüdeki gibi),
yakınsama ayrıca tutarlı bir dren rayı ister. Testler: `tests/test_exp20.py`; tüm takım 74/74.

**Örtük çözücüyle doğrulama (`results/exp20/validate`):** Exp19 basamakları 18/20 aynı k; halka 19/20 aynı; ani öldürme
doğumları 3,3 V'ta 2 → 2, 5,0 V'ta 2 → 3; Exp17 taramalarında ilk bölünme 5 koşunun 4'ünde aynı (E_lam60: 3,34 V/101 µm →
3,42 V/115,5 µm, öngörü 108 µm); sonraki bölünmeler eşik yakınında yola bağlı (%15'e kadar fark).

**Örtük çözücüyle çukurlar (54 koşu, `predictions_dip_imp.json`):**

| tahmin | sonuç |
|---|---|
| D1: V_DD çukuru, derinlik hücre payını (3,3 V: 7,8–10,2 mV; 5,0 V: 2,7–3,2 mV) aşınca böler | ✗ — 0,1 V'a kadar (payın 10–37 katı) hiç bölünme yok; 5,0 V'ta 0,2 V: 3 → 4; 0,5 V: 3,3 V'ta 3 → 5, 5,0 V'ta 3 → 6. Süre (0,02–0,5 s) fark etmiyor. |
| D2: kapı çukuru sırasında bölünme yok | ✓ — 12 kapı koşusunun hiçbirinde çukur içinde bölünme yok |
| D3: 0,5 V çukurlar öldürür, bölmez | ✗ — böldürüyor, ölüm yok |
| D4: kapı çukurunun bırakılması böler, derinlikle artar, 5,0 V'ta daha çok; 5,0 V −60 mV: iki sürede de 3 → ≥ 6 | ✓ — bölünmeler bırakmadan 20 ms sonra; 3,3 V: −20/−40/−60 mV → 3/5/6 hücre; 5,0 V → 5/6/7; −60 mV'ta 3 → 7 (iki süre) |

Sonradan (kayıtsız) ortak tablo: toplam akımın ani zorla küçülmesi (çukurun bırakılması ya da ≥ 0,2 V V_DD düşüşünün ilk
20 ms'i) cepheler geri çekilemeden merkezleri V_down altına itiyor; hücreler eşzamanlı bölünüyor. D1'in eşiği yanlıştı çünkü
"pay" sabit j_on varsayıyordu; gerilim düşünce tüm hücrelerin akımı birlikte azalıyor ve merkez–kenar farkı da küçülüyor.

## Exp17 yeniden sınaması (örtük çözücü, 40 koşu, `recheck17/analysis.json`)

| büyüklük | yayımlanan | örtük çözücü |
|---|---|---|
| mekanizma oranı, temiz | 1,001 ± 0,069 (n = 22) | 1,029 ± 0,071 (n = 20) |
| mekanizma oranı, sert yük | 1,038 ± 0,062 (n = 22) | 1,003 ± 0,057 (n = 24) |
| mekanizma oranı, evrensellik | 0,979 ± 0,133 (n = 9) | 0,973 ± 0,116 (n = 9) |
| ilk bölünme, temiz (N seti hariç) | n = 14, ort. +%1,6, RMS %9,1 | n = 13, ort. +%3,2, RMS %8,6 |

- 14 yayımlanmış ilk bölünmenin 13'ü yeniden üretiliyor. G_rD5_lam80 (eskiden 5,38 V'ta %+20) artık tarama sonuna (5,6 V)
  kadar bölünmüyor.
- Geçerlilik dışı düşük sızıntı koşusu (N_ballast_lam60, g_sub = 1e-15): yayımlanan metin "1 s/adımda bölünme yok" diyor;
  örtük çözücüde 175 µm'de bölünüyor (öngörünün ~2 katı, mekanizma oranı ~2,1). "Geçerlilik dışı" sınıflaması değişmiyor,
  o cümle değişmeli.
- Exp17 raporunun başına bu bulguyu anlatan "Düzeltme 2" kutusu eklendi.

## Yapılmayan

- E (hücre döngüsü 1 → 2 → 4 → 8): zaman kalmadığı için yapılmadı. D4 ve V_DD çukuru sonuçları buna doğrudan yol veriyor.
