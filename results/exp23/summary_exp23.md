# Exp23 — Programlanabilir hücrenin geri dönüştürülen enerjisindeki kayıp ve hücrelerin atom kontrolüne etkisi

Tarih: 2026-10-15. Kod: `experiments/exp23_energy_atoms/e23.py` (`python e23.py 3.3`, `python e23.py 5`, `fig`).
Cihaz: Exp18 durumları (3 hücre, tohum 1), örtük çözücü. Şekil: `fig_exp23.png`.

## 1. Geri dönüştürülen enerjideki kayıp

**Yöntem:**
- En sağdaki hücre (−300 mV kapı bandı, 30 µm kenar payı) iki şekilde öldürüldü: yavaş (30 s rampa) ve ani (20 ms adımlar).
- Tam profillerden güç ayrıştırıldı: kaynak = yük direnci + balastlar + dren rayı + kanallar (hücreler).
- Muhasebe kapanıyor: hata 0,000 µW.

| durum | hücre | ölen hücrenin gücü | kalanlara geçen | verim η | geçişin ek enerjisi |
|---|---|---|---|---|---|
| 3,3 V yavaş | 3 → 2 | 135,9 µW | 133,9 µW | 0,985 | 16,3 µJ (35 s) |
| 3,3 V ani | 3 → 4 | 135,9 µW | 137,2 µW | 1,009 | ≈ 0 |
| 5,0 V yavaş | 3 → 3 (kendiliğinden yenilenme) | 220,0 µW | 222,2 µW | 1,010 | 93,5 µJ (35 s) |
| 5,0 V ani | 3 → 5 | 220,0 µW | 225,0 µW | 1,023 | 0,35 µJ |

- **Geri dönüşüm neredeyse kayıpsız.** Toplam akımı yük doğrusu belirliyor; ölen hücrenin akımı kalanlara geçmek zorunda.
- **En kötü durum (3,3 V yavaş):** %1,5 kayıp. Hücrelerin toplam gücü 392,9 → 390,2 µW (−%0,7); balast +3,7 µW,
  ray +0,9 µW ısınıyor.
- **Ani ölüm yeni hücre doğuruyor.** Bu, Exp18–19 ile tutarlı.
- **Asıl kayıp geri dönüşümde değil, mimaride.** Kaynak gücünün %74,6'sı (3,3 V) ve %83'ü (5 V) her durumda yük
  direncinde yanıyor. Hücrelere giden pay %24 ve %16.
- **Ölen hücrenin gövde enerjisi** (½ c_b V_b²) yalnızca 3–5e-14 J. Geri kazanılacak depolanmış enerji yok; hücre bir
  "akış" yapısı, enerji deposu değil.

## 2. Atomları programlayabilir miyiz? (test)

Simüle edilen durumların büyüklükleri, atomu yerinden oynatma ya da değiştirme eşikleriyle karşılaştırıldı:

| büyüklük | hücrelerde | gereken | oran |
|---|---|---|---|
| taşıyıcı enerjisi (tepe qV_D) | 0,82 eV (3,3 V ve 5 V'ta aynı) | Si–H bağı 3,6 eV; kafes atomu sökme 13 eV | 0,23; 0,06 |
| e–e saçılma kuyruğu (2qV_D) | 1,64 eV | Si–H 3,6 eV | 0,46 |
| yanal alan | 0,02–0,07 V/nm | STM ile atom itme ~1 V/nm; alan buharlaşması ~10 V/nm | ≤ 0,07 |
| dikey oksit alanı | ~0,12 V/nm | iyon sürüklenmesi 4 K'de donmuş | — |
| ısınma (çizgi ısıtıcı) | 0,03–5 K | katkı difüzyonu ~600 °C | ≤ 0,01 |
| akımın manyetik alanı (1 µm'de) | 24 µT, Zeeman 1,4e-9 eV | kT(4 K) = 3,4e-4 eV | 4e-6 |
| sıcak taşıyıcı hasar hızı | nominal 1,8 V'un 3e-11 katı | — | pratikte sıfır |
| kanal/kontak/li1 akım yoğunluğu | ~1e5 A/cm² | elektromigrasyon (300 K) ~1e5 A/cm²; 4 K'de ısıl olarak donmuş | — |
| uzaysal çözünürlük | hücre ~40 µm | atom ~0,5 nm | 1e-5 |

- **Keşif: mimarinin yapısal bir özelliği.** Filaman içindeki gerilim bistabilite tarafından V_M civarına kenetleniyor.
  V_DD 3,3 → 5 V'a çıkarılsa da taşıyıcı enerjisi artmıyor (0,82 eV); fazla güç yalnızca yeni hücre doğuruyor. Bu bir
  başarısızlık değil, mimarinin doğasının bulunması: hücre enerjiyi yoğunlaştırmaz, sayıya çevirir. Örgü atomunu kendi
  enerjisiyle oynatamamasının nedeni de budur.
- **Olumlu yetenek: toplu yük durumu kontrolü.** Hücrenin açılıp kapanması gövdede µm başına ~4400 bor alıcısının yük durumunu
  (B⁻ ↔ B⁰) topluca, geri dönüşlü ve µm konum seçimiyle değiştiriyor (c_b × 0,7 V / q). Bu klasik bir işlem için yeterli.
  Tek başına "atom programlama" değil. Exp24 bunu tek atom düzeyine indirdi: aynı gövdede atomların her biri kendi
  gerilim adresinde tek tek 0/1 yapılabiliyor (model düzeyinde). 4 K'de bunun için ~60 nm gövde ve dengede kalibre adres
  gerekiyor (~%99 doğruluk, gövde başına ~3 atom); 1 K'de 100 nm gövdede ~35 atom.
- **Tek atom adreslemesi:** en küçük SKY130 kanalında ~750 katkı atomu var (atom başına ~5900 nm³). Tek atom için
  ~10–20 nm'lik nanotel ve ≤ 4 K gerekir.
- **300 K'de yüksek akım:** atomlara yalnızca rastgele ve yıkıcı etki eder (temaslarda elektromigrasyon boşlukları,
  arayüz tuzakları). Bu "programlama" değil, aşınma.

**Sonuç (Exp24 ile güncellendi):**
- Hücrenin kendi enerjisiyle örgü atomunu doğrudan oynatmak mümkün değil; tüm eşiklerin 2 ile 11 mertebe altındayız.
- Hücre gövdesindeki doğal atomların yük durumu topluca kontrol ediliyor (bu deney).
- Aynı atomlar gerilim adresiyle tek tek 0/1 yapılabiliyor; atom konumu ve spini ise dolaylı yollarla programlanabiliyor
  (ReRAM zinciri, tek atomluk temas, ³¹P spini). Ayrıntı: `results/exp24/summary_exp24.md`.

## Sınırlar

- Enerji muhasebesi indirgenmiş 4 K modelinde (Exp18 durumları); geçiş enerjisi adım sonu örneklerinden.
- Atom tablosu mertebe tahmini. Eşik değerleri literatürdeki tipik değerler (Si–H 3,6 eV, STM ~1 V/nm, elektromigrasyon
  ~1e5 A/cm²); burada ayrıca doğrulanmadı.
