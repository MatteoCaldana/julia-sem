#pragma once

namespace sem {
  
double leg___00(double x0) {
  return 1.;
}

double leg_p_00(double x0) {
  return 0.;
}

double legpp_00(double x0) {
  return 1.;
}

double leg___01(double x0) {
  return x0;
}

double leg_p_01(double x0) {
  return 1.;
}

double legpp_01(double x0) {
  return 0.;
}

double leg___02(double x0) {
  const auto x1 = x0 * x0;
  return (-1 + (3.) * x1 ) / 2.;
}

double leg_p_02(double x0) {
  return 3. * x0 ;
}

double legpp_02(double x0) {
  return 3.;
}

double leg___03(double x0) {
  const auto x1 = x0 * x0;
  return ( (-3.) * x0 + (5.) * x0*x1 ) / 2.;
}

double leg_p_03(double x0) {
  const auto x1 = x0 * x0;
  return (-3 + (15.) * x1 ) / 2.;
}

double legpp_03(double x0) {
  return 15. * x0 ;
}

double leg___04(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (3 + (-30.) * x1 + (35.) * x2 ) / 8.;
}

double leg_p_04(double x0) {
  const auto x1 = x0 * x0;
  return ( (-15.) * x0 + (35.) * x0*x1 ) / 2.;
}

double legpp_04(double x0) {
  const auto x1 = x0 * x0;
  return (-15 + (105.) * x1 ) / 2.;
}

double leg___05(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (15.) * x0 + (-70.) * x0*x1 + (63.) * x0*x2 ) / 8.;
}

double leg_p_05(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (15 + (-210.) * x1 + (315.) * x2 ) / 8.;
}

double legpp_05(double x0) {
  const auto x1 = x0 * x0;
  return ( (-105.) * x0 + (315.) * x0*x1 ) / 2.;
}

double leg___06(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (-5 + (105.) * x1 + (-315.) * x2 + (231.) * x1*x2 ) / 16.;
}

double leg_p_06(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (105.) * x0 + (-630.) * x0*x1 + (693.) * x0*x2 ) / 8.;
}

double legpp_06(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (105 + (-1890.) * x1 + (3465.) * x2 ) / 8.;
}

double leg___07(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (-35.) * x0 + (315.) * x0*x1 + (-693.) * x0*x2 + (429.) * x0*x1*x2 ) / 16.;
}

double leg_p_07(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (-35 + (945.) * x1 + (-3465.) * x2 + (3003.) * x1*x2 ) / 16.;
}

double legpp_07(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (945.) * x0 + (-6930.) * x0*x1 + (9009.) * x0*x2 ) / 8.;
}

double leg___08(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (35 + (-1260.) * x1 + (6930.) * x2 + (-12012.) * x1*x2 + (6435.) * x3 ) / 128.;
}

double leg_p_08(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (-315.) * x0 + (3465.) * x0*x1 + (-9009.) * x0*x2 + (6435.) * x0*x1*x2 ) / 16.;
}

double legpp_08(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return (-315 + (10395.) * x1 + (-45045.) * x2 + (45045.) * x1*x2 ) / 16.;
}

double leg___09(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (315.) * x0 + (-4620.) * x0*x1 + (18018.) * x0*x2 + (-25740.) * x0*x1*x2 + (12155.) * x0*x3 ) / 128.;
}

double leg_p_09(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (315 + (-13860.) * x1 + (90090.) * x2 + (-180180.) * x1*x2 + (109395.) * x3 ) / 128.;
}

double legpp_09(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  return ( (-3465.) * x0 + (45045.) * x0*x1 + (-135135.) * x0*x2 + (109395.) * x0*x1*x2 ) / 16.;
}

double leg___10(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-63 + (3465.) * x1 + (-30030.) * x2 + (90090.) * x1*x2 + (-109395.) * x3 + (46189.) * x1*x3 ) / 256.;
}

double leg_p_10(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (3465.) * x0 + (-60060.) * x0*x1 + (270270.) * x0*x2 + (-437580.) * x0*x1*x2 + (230945.) * x0*x3 ) / 128.;
}

double legpp_10(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (3465 + (-180180.) * x1 + (1351350.) * x2 + (-3063060.) * x1*x2 + (2078505.) * x3 ) / 128.;
}

double leg___11(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-693.) * x0 + (15015.) * x0*x1 + (-90090.) * x0*x2 + (218790.) * x0*x1*x2 + (-230945.) * x0*x3 + (88179.) * x0*x1*x3 ) / 256.;
}

double leg_p_11(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-693 + (45045.) * x1 + (-450450.) * x2 + (1531530.) * x1*x2 + (-2078505.) * x3 + (969969.) * x1*x3 ) / 256.;
}

double legpp_11(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (45045.) * x0 + (-900900.) * x0*x1 + (4594590.) * x0*x2 + (-8314020.) * x0*x1*x2 + (4849845.) * x0*x3 ) / 128.;
}

double leg___12(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (231 + (-18018.) * x1 + (225225.) * x2 + (-1021020.) * x1*x2 + (2078505.) * x3 + (-1939938.) * x1*x3 + (676039.) * x2*x3 ) / 1024.;
}

double leg_p_12(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-9009.) * x0 + (225225.) * x0*x1 + (-1531530.) * x0*x2 + (4157010.) * x0*x1*x2 + (-4849845.) * x0*x3 + (2028117.) * x0*x1*x3 ) / 256.;
}

double legpp_12(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-9009 + (675675.) * x1 + (-7657650.) * x2 + (29099070.) * x1*x2 + (-43648605.) * x3 + (22309287.) * x1*x3 ) / 256.;
}

double leg___13(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (3003.) * x0 + (-90090.) * x0*x1 + (765765.) * x0*x2 + (-2771340.) * x0*x1*x2 + (4849845.) * x0*x3 + (-4056234.) * x0*x1*x3 + (1300075.) * x0*x2*x3 ) / 1024.;
}

double leg_p_13(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (3003 + (-270270.) * x1 + (3828825.) * x2 + (-19399380.) * x1*x2 + (43648605.) * x3 + (-44618574.) * x1*x3 + (16900975.) * x2*x3 ) / 1024.;
}

double legpp_13(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-135135.) * x0 + (3828825.) * x0*x1 + (-29099070.) * x0*x2 + (87297210.) * x0*x1*x2 + (-111546435.) * x0*x3 + (50702925.) * x0*x1*x3 ) / 256.;
}

double leg___14(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-429 + (45045.) * x1 + (-765765.) * x2 + (4849845.) * x1*x2 + (-14549535.) * x3 + (22309287.) * x1*x3 + (-16900975.) * x2*x3 + (5014575.) * x1*x2*x3 ) / 2048.;
}

double leg_p_14(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (45045.) * x0 + (-1531530.) * x0*x1 + (14549535.) * x0*x2 + (-58198140.) * x0*x1*x2 + (111546435.) * x0*x3 + (-101405850.) * x0*x1*x3 + (35102025.) * x0*x2*x3 ) / 1024.;
}

double legpp_14(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (45045 + (-4594590.) * x1 + (72747675.) * x2 + (-407386980.) * x1*x2 + (1003917915.) * x3 + (-1115464350.) * x1*x3 + (456326325.) * x2*x3 ) / 1024.;
}

double leg___15(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-6435.) * x0 + (255255.) * x0*x1 + (-2909907.) * x0*x2 + (14549535.) * x0*x1*x2 + (-37182145.) * x0*x3 + (50702925.) * x0*x1*x3 + (-35102025.) * x0*x2*x3 + (9694845.) * x0*x1*x2*x3 ) / 2048.;
}

double leg_p_15(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-6435 + (765765.) * x1 + (-14549535.) * x2 + (101846745.) * x1*x2 + (-334639305.) * x3 + (557732175.) * x1*x3 + (-456326325.) * x2*x3 + (145422675.) * x1*x2*x3 ) / 2048.;
}

double legpp_15(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (765765.) * x0 + (-29099070.) * x0*x1 + (305540235.) * x0*x2 + (-1338557220.) * x0*x1*x2 + (2788660875.) * x0*x3 + (-2737957950.) * x0*x1*x3 + (1017958725.) * x0*x2*x3 ) / 1024.;
}

double leg___16(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return (6435 + (-875160.) * x1 + (19399380.) * x2 + (-162954792.) * x1*x2 + (669278610.) * x3 + (-1487285800.) * x1*x3 + (1825305300.) * x2*x3 + (-1163381400.) * x1*x2*x3 + (300540195.) * x4 ) / 32768.;
}

double leg_p_16(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-109395.) * x0 + (4849845.) * x0*x1 + (-61108047.) * x0*x2 + (334639305.) * x0*x1*x2 + (-929553625.) * x0*x3 + (1368978975.) * x0*x1*x3 + (-1017958725.) * x0*x2*x3 + (300540195.) * x0*x1*x2*x3 ) / 2048.;
}

double legpp_16(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return (-109395 + (14549535.) * x1 + (-305540235.) * x2 + (2342475135.) * x1*x2 + (-8365982625.) * x3 + (15058768725.) * x1*x3 + (-13233463425.) * x2*x3 + (4508102925.) * x1*x2*x3 ) / 2048.;
}

double leg___17(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return ( (109395.) * x0 + (-5542680.) * x0*x1 + (81477396.) * x0*x2 + (-535422888.) * x0*x1*x2 + (1859107250.) * x0*x3 + (-3650610600.) * x0*x1*x3 + (4071834900.) * x0*x2*x3 + (-2404321560.) * x0*x1*x2*x3 + (583401555.) * x0*x4 ) / 32768.;
}

double leg_p_17(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return (109395 + (-16628040.) * x1 + (407386980.) * x2 + (-3747960216.) * x1*x2 + (16731965250.) * x3 + (-40156716600.) * x1*x3 + (52933853700.) * x2*x3 + (-36064823400.) * x1*x2*x3 + (9917826435.) * x4 ) / 32768.;
}

double legpp_17(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  return ( (-2078505.) * x0 + (101846745.) * x0*x1 + (-1405485081.) * x0*x2 + (8365982625.) * x0*x1*x2 + (-25097947875.) * x0*x3 + (39700390275.) * x0*x1*x3 + (-31556720475.) * x0*x2*x3 + (9917826435.) * x0*x1*x2*x3 ) / 2048.;
}

double leg___18(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return (-12155 + (2078505.) * x1 + (-58198140.) * x2 + (624660036.) * x1*x2 + (-3346393050.) * x3 + (10039179150.) * x1*x3 + (-17644617900.) * x2*x3 + (18032411700.) * x1*x2*x3 + (-9917826435.) * x4 + (2268783825.) * x1*x4 ) / 65536.;
}

double leg_p_18(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return ( (2078505.) * x0 + (-116396280.) * x0*x1 + (1873980108.) * x0*x2 + (-13385572200.) * x0*x1*x2 + (50195895750.) * x0*x3 + (-105867707400.) * x0*x1*x3 + (126226881900.) * x0*x2*x3 + (-79342611480.) * x0*x1*x2*x3 + (20419054425.) * x0*x4 ) / 32768.;
}

double legpp_18(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return (2078505 + (-349188840.) * x1 + (9369900540.) * x2 + (-93699005400.) * x1*x2 + (451763061750.) * x3 + (-1164544781400.) * x1*x3 + (1640949464700.) * x2*x3 + (-1190139172200.) * x1*x2*x3 + (347123925225.) * x4 ) / 32768.;
}

double leg___19(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return ( (-230945.) * x0 + (14549535.) * x0*x1 + (-267711444.) * x0*x2 + (2230928700.) * x0*x1*x2 + (-10039179150.) * x0*x3 + (26466926850.) * x0*x1*x3 + (-42075627300.) * x0*x2*x3 + (39671305740.) * x0*x1*x2*x3 + (-20419054425.) * x0*x4 + (4418157975.) * x0*x1*x4 ) / 65536.;
}

double leg_p_19(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return (-230945 + (43648605.) * x1 + (-1338557220.) * x2 + (15616500900.) * x1*x2 + (-90352612350.) * x3 + (291136195350.) * x1*x3 + (-546983154900.) * x2*x3 + (595069586100.) * x1*x2*x3 + (-347123925225.) * x4 + (83945001525.) * x1*x4 ) / 65536.;
}

double legpp_19(double x0) {
  const auto x1 = x0 * x0;
  const auto x2 = x1 * x1;
  const auto x3 = x2 * x2;
  const auto x4 = x3 * x3;
  return ( (43648605.) * x0 + (-2677114440.) * x0*x1 + (46849502700.) * x0*x2 + (-361410449400.) * x0*x1*x2 + (1455680976750.) * x0*x3 + (-3281898929400.) * x0*x1*x3 + (4165487102700.) * x0*x2*x3 + (-2776991401800.) * x0*x1*x2*x3 + (755505013725.) * x0*x4 ) / 32768.;
}


constexpr double(*legendre[])(double) = {
  &leg___00,
  &leg___01,
  &leg___02,
  &leg___03,
  &leg___04,
  &leg___05,
  &leg___06,
  &leg___07,
  &leg___08,
  &leg___09,
  &leg___10,
  &leg___11,
  &leg___12,
  &leg___13,
  &leg___14,
  &leg___15,
  &leg___16,
  &leg___17,
  &leg___18,
  &leg___19
};

constexpr double(*legendre_prime[])(double) = {
  &leg_p_00,
  &leg_p_01,
  &leg_p_02,
  &leg_p_03,
  &leg_p_04,
  &leg_p_05,
  &leg_p_06,
  &leg_p_07,
  &leg_p_08,
  &leg_p_09,
  &leg_p_10,
  &leg_p_11,
  &leg_p_12,
  &leg_p_13,
  &leg_p_14,
  &leg_p_15,
  &leg_p_16,
  &leg_p_17,
  &leg_p_18,
  &leg_p_19
};

constexpr double(*legendre_second[])(double) = {
  &legpp_00,
  &legpp_01,
  &legpp_02,
  &legpp_03,
  &legpp_04,
  &legpp_05,
  &legpp_06,
  &legpp_07,
  &legpp_08,
  &legpp_09,
  &legpp_10,
  &legpp_11,
  &legpp_12,
  &legpp_13,
  &legpp_14,
  &legpp_15,
  &legpp_16,
  &legpp_17,
  &legpp_18,
  &legpp_19
};
}