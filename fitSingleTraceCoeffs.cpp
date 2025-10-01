#include <vector>
#include "FAWSolver.h" // Assuming this is the header for FAWSolver
#include "WaveformAnalyzer.h" // Assuming this is the header for WaveformAnalyzer
#include "SeisTraceParams.h" // Assuming this is the header for SeisTraceParams
#include <QtCore/qmath.h>
#include <QMap>
#include <Eigen/Dense>

/*
 * Note: This is a member function of the ProcNextgenFirstBreakPicker class.
 * It relies on member variables such as m_dt, m_tuneToEnergyCenter,
 * m_signalWidthCoefficent, m_polynomDegree, m_lambda, m_weightingMode,
 * m_degreePenalty, m_compressionCoefficient, coeffsMap, and traceParamsMap.
 *
 * This file is provided to show the standalone implementation of the function
 * as requested.
 */
void ProcNextgenFirstBreakPicker::fitSingleTraceCoeffs(int traceId,
							  WaveformAnalyzer& wfAnal,
							  int lengthOfBuffer,
							  float rmsSignalDuration,
							  int startIndex,
							  std::vector<float>& signalToAutocorr,
							  std::vector<float>& autocorrRes,
							  std::vector<float>& energyCenters)
{
    float* bufferPtr = wfAnal.getProcessedBufferForTrace(traceId, startIndex);
    doFrontSmooth(bufferPtr, lengthOfBuffer);
    doBackSmooth(bufferPtr, lengthOfBuffer);
    int tempEnergyCenter = qFloor(findApproximateEnergyCenter(bufferPtr, 0, lengthOfBuffer));

    float energyCenterTime;
    if (startIndex < 0) {
        energyCenterTime = static_cast<float>(qAbs(tempEnergyCenter + startIndex) * m_dt);
    } else {
        energyCenterTime = static_cast<float>((tempEnergyCenter + startIndex) * m_dt);
    }
    energyCenters.push_back(energyCenterTime);

    memcpy(&signalToAutocorr.at(0), bufferPtr, sizeof(float) * lengthOfBuffer);
    autocorrRes = calculateAutocorrelation_FFT(signalToAutocorr);
    int tempSignalDuration = 0;

    if (m_tuneToEnergyCenter) {
        tempSignalDuration = qFloor(rmsSignalDuration / 2);
    } else {
        tempSignalDuration = findDominantPeriod(autocorrRes, 1);
    }

    if (m_tuneToEnergyCenter) {
        doFrontSmooth(bufferPtr, tempEnergyCenter, tempSignalDuration);
    } else {
        doFrontSmooth(bufferPtr, tempEnergyCenter, tempSignalDuration * m_signalWidthCoefficent);
    }

    Eigen::VectorXf tempCoefs;
    if (m_tuneToEnergyCenter) {
        tempCoefs = FAWSolver::fitLegendreApproximationExtended(bufferPtr, lengthOfBuffer, qFloor(tempEnergyCenter),
                                                                qFloor(tempSignalDuration),
                                                                m_polynomDegree, m_lambda, m_weightingMode,
                                                                m_degreePenalty, m_compressionCoefficient);
    } else {
        tempCoefs = FAWSolver::fitLegendreApproximationExtended(bufferPtr, lengthOfBuffer, qFloor(tempEnergyCenter),
                                                                qFloor(tempSignalDuration) * m_signalWidthCoefficent,
                                                                m_polynomDegree, m_lambda, m_weightingMode,
                                                                m_degreePenalty, m_compressionCoefficient);
    }
    coeffsMap.insert(traceId, tempCoefs);
    traceParamsMap[traceId].windowCenterTime = energyCenterTime;
    if (m_tuneToEnergyCenter) {
        traceParamsMap[traceId].windowHalfWidth = qFloor(tempSignalDuration);
    } else {
        traceParamsMap[traceId].windowHalfWidth = qFloor(tempSignalDuration) * m_signalWidthCoefficent;
    }
}