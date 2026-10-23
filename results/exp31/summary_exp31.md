# Exp31 — Hücre yasalarının ölçülmüş akım filamanlarıyla sınanması ve bütün verinin son kapanışı

Tarih: 2026-10-16. Kod: `experiments/exp31_cell_laws_measured/cl31.py` (`python cl31.py fig`). Öngörüler hesaptan önce
yazıldı: `predictions_exp31.json` (15:14). Çıktılar:
- `cell_laws.json`;
- `closure_final.json` (Exp30'un 28 satırı, Exp31 güncellemeleriyle);
- `fig_exp31.png`.

**Neden.** Stop kancası: Exp30'dan sonra ölçülmemiş kalan tek katman silisyum hücre yasalarıydı (Exp13, 14, 17, 19, 20, 23). Laboratuvar
yerine, aynı fizik sınıfının (S biçimli negatif direnç + küresel bağlaşım + balast) silisyumda ölçülmüş yayımlanmış verisi kullanıldı.

**Açıklama (kör değil).** Tezi okurken bazı oranları kafadan kabaca hesapladım. Sınanan yasalar Exp17'de çok önce sabitlenmişti ve
bu veriye göre ayarlanmadı. Yöntem önceden sabitlendi; olumsuz sonuç (CL2) dahil hepsi yazıldı.

**Fizik farkı.** Ölçülen filamanlar 300 K'de, bipolar geri dönüşlü (snapback), ısıl kökenli ve darbeli (TLP, ~100 ns–1 µs).
Projenin hücreleri 4 K'de, darbe iyonlaşması ve donmuş gövdeyle. Sınanan şey yasaların yapısı; SKY130 sabitleri değil.

## Kısa sonuç

- **Exp17 hücre yasası ölçülmüş silisyum filamanlarında büyük ölçüde tutuyor** (Mamanee 2012, TU Wien, ns/µm çözünürlüklü girişim
  haritalaması):
  - sabit filaman sayısında filaman başına akım yoğunluğu sabit, genişlik akımla büyüyor (%0,1 ve %5,2);
  - sayı yasası N = I/(0,75 j w_s) ölçülen aralıklara düşüyor (ölçülen ortalamanın 0,84–0,90'ı);
  - tek filamanın en büyük akımı j·w_s = 0,889 A, ölçülen bölünme penceresinin üst kenarı 0,9 A (−%1,2);
  - **tutmayan:** yoğunluk filaman sayısıyla %22 düşüyor (Exp17'de tek bir j_on vardı).
- **Programlama, bölünme ve paylaşım ölçülmüş:**
  - tohumsuz aygıtta doğum yeri rastgele, yerel katkı tohumlu aygıtta sabit (Exp18);
  - filaman ikiye bölünüyor ve güç eşit paylaşılıyor (Exp13, 23);
  - kısa yüksek darbede daha çok filaman çıkıyor (Exp19);
  - V₂O₃'te akım düşürülünce filaman çok sayıda paralel filamana bölünüyor (Exp20; mekanizma farklı);
  - aktif damlalar içi boş sıvı kabuğa dönüşüyor (Exp14; bölünme penceresi ölçülmedi).
- **Son sayım (Exp0–27, 28 sonuç):** 13 ✓, 7 karışık/kısmen, 5 ✗, 3 yöntem, **0 ölçüm yok**.

## Nicel sınamalar (Mamanee 2012, Tablo 8 ve Şekil 73–74)

Ölçülen filaman genişlikleri W_CF ve bunlardan j = I/(N·W):

| I_S | N | W_CF (µm) | j (mA/µm) |
|---|---|---|---|
| 1,4 A | 2 | 15–16,6 | 42,2–46,7 |
| 1,4 A | 3 | 12,7–13,4 | 34,8–36,7 |
| 2 A | 2 | 22–23 | 43,5–45,5 |
| 2 A | 3 | 16,9–18,6 | 35,8–39,4 |
| 2 A | 4 | 13,6–15,2 | 32,9–36,8 |

| no | yasa | sonuç |
|---|---|---|
| CL1 | sabit N'de j sabit (Exp17: w = I/(N j_on)) | ✓: N = 2 için +%0,1; N = 3 için +%5,2 (ölçüt ≤ %10) |
| CL2 | j, N'den bağımsız | **✗**: N = 2 → 4'te −%21,7 |
| CL3 | N = I/(0,75 j w_s), j_ref = 44,4 mA/µm, w_s = 20 µm | ✓: 0,9 / 1,4 / 2 A'de 1,35 / 2,10 / 3,00; ölçülen 1–2 / 2–3 / 3–4. w_s = 22,5 µm ile de aralıkta |
| CL4 | tek filaman sınırı j·w_s | ✓: 0,889 A vs pencerenin üst kenarı 0,9 A. Alt kenar (0,45 A'den itibaren iki filaman mümkün) Exp17'de öngörülmüyor: çok-kararlılık |
| CL5 | bölünme ve eşit güç paylaşımı (Exp13, 23) | ✓: 0,8 A'de filaman ikiye bölünüyor, iki filaman gücü eşit paylaşıyor |
| CL6 | tohum doğum yerini seçer (Exp18) | ✓: tohumsuz aygıtlarda I_S < 0,45 A'de rastgele başlangıç; NPN_Dope'ta (P-gövde ortasında yüksek katkı) ilk filaman tekrarlanabilir biçimde ortada |
| CL7 | ani/kısa sürüş daha çok hücre (Exp19) | ✓: kısa darbede yüksek akımda daha çok filaman; sabit akımda N darbeden darbeye değişiyor |
| CL8 | gerilim kenedi ve balast (Exp12) | ✓ (yapı): I–V'de her N için ayrı dal; iki filamanda tutma gerilimi ~1 V düşük, ΔV = I R_A (1 − 1/N). SKY130'un V_M = 0,80 V'u ölçülmedi |
| CL9 | 300 K'de kalan yol bipolar (Exp21) | ✓: 300 K'de bipolar filamanlar çoğalıyor ve bölünüyor (darbeli; yüksek akımda ısıl bozulma) |
| CL10 | akım düşüşüyle bölünme (Exp20) | kısmen: V₂O₃'te akım düşürülünce filaman önce dallanıyor, sonra çok sayıda paralel filamana bölünüyor, sonra filamanlar tek tek sönüyor. Yazarlar bunu gerinim heterojenliğine bağlıyor (Lange 2020) |
| CL11 | boyut penceresi (Exp14) | kısmen: aktif koaservat damlaları kritik hacmin üstünde içi boş sıvı kabuğa dönüşüyor (Bergmann 2023), Exp14'ün pencere üstü öngörüsü. Aktif damla bölünmesi henüz deneyle gösterilmemiş |

## Bütün verinin son kapanışı (Exp30 + Exp31)

| sonuç | sayı | deneyler |
|---|---|---|
| ✓ | 13 | Exp3, 4, 5, 8, 11, 12, 13, 15, 18, 19, 21, 23, 27 |
| karışık / kısmen | 7 | Exp0, 2, 10, 14, 17, 20, 24 |
| ✗ | 5 | Exp1 (ada kapasitansı 1,7–3,8×), Exp6 (eşleşmezlik 3,5–5,6×), Exp22 (alt bant), Exp25 (GHZ), Exp26 (eşevrelilik) |
| yöntem | 3 | Exp7, 9, 16 |
| ölçüm yok | 0 | — |

Satır satır tablo: `closure_final.json`.

## Doğru yol ve denklem (ölçülmüş dayanaklarıyla)

  F = γ · I_toplam · √(2 m V / q),   I_toplam = Σ hücre akımları = N(I) · I_c(N)

| terim | ölçülmüş dayanak |
|---|---|
| F = γ I √(2mV/q) | NEXT 40 nokta, kök ortalama kare %1,45 (Exp29); In FEEP, doğrudan itki %1,6–4,3 (Exp30) |
| alanla açılan atom yayımı, yayıcı sayısının akımla artması | FEEP iğneleri 1/28 → 28/28; SiOx programlanabilir yayıcı (Exp29) |
| tek yük / tek atom metrelemesi | Si pompası 0,27 ppm; istek üzerine tek iyon 123/139 (Exp29) |
| hücreler: sabit j'de genişleme, sayının akımla basamaklı artması, bölünme, tohumla doğum, eşit paylaşım | silisyum filamanları (Mamanee 2012); V₂O₃ (Lange 2020) |
| kriyojenik sınıf sınırı | 180 nm CMOS kink 30 K'de var, 40 K'de yok (Exp30) |
| turnike adası | ölçümde ~2 nm'de 300 K (Exp30) |

**Önemli sonuç.** İtki yalnız toplam akıma bağlı ve bu bağ ölçümle doğrulandı. Hücre yasasındaki kusur (CL2) itki öngörüsünü değiştirmiyor;
yalnız toplam akımın kaç hücreye bölündüğünü (basamak yapısını) değiştiriyor. Düzeltilmiş yasa:

  N ≈ I / (0,75 · j(N) · w_s),   j(N) N ile azalır (ölçülen: N = 2 → 4'te −%22)

## Sınırlar

- **SKY130'a özgü sabitler ölçülmedi:** 4 K'de j_on = 2,569 µA/µm ve V_M = 0,80 V. Ölçülen filamanlar başka bir aygıt sınıfında
  (300 K bipolar, darbeli). Gerekli ölçüm: ≤ 30 K'de geniş SKY130 nFET ve bölümlenmiş drenaj (Exp12/17 önerileri).
- **Darbeli rejim.** Ölçülen bipolar filamanlar darbeli rejimde ve yüksek akımda ısıl bozulmaya gidiyor. Sürekli akım kaynağı olarak
  ısıl yönetim gerekir.
- **Exp20.** Olgu ölçülmüş, mekanizma farklı.
- **Exp14.** Bölünme penceresinin kendisi ölçülmedi.

## Kaynaklar

- W. Mamanee, *Investigation of current filament behavior in smart power ESD devices*, TU Wien tezi (2012, danışman D. Pogany),
  https://repositum.tuwien.at/handle/20.500.12708/13309 (tam metin; §3.2–3.3, Tablo 8)
- M. Lange vd., *Optical imaging of strain-mediated phase coexistence during electrothermal switching in a Mott insulator*,
  arXiv 2009.12536 (tam metin)
- A. M. Bergmann vd., *Liquid spherical shells are a non-equilibrium steady state of active droplets*, Nat. Commun. (2023),
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10582082/ (özet)
