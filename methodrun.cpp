void ProcNextgenFirstBreakPicker::run()
{
	DatasetPassport* passp = getPassport();
	QStringList missedWords = checkNecessaryWords(*passp);
	if(!missedWords.isEmpty()){
		QString message = QString("Unable to create FirstBreaks object. Necessary words are missing (%1).").arg(missedWords.join(", "));
		send( message,MessageEvent::Error );
		return ;
	}
	TraceHeaderMap* thm = getTraceHeaderMap();

	m_nSamples = passp->getNumSamples();
	m_dt = passp->getSampleInterval() / 1000.0f;

	if (!guidPrep(passp, thm)){
		send(QString("Cannot read input file"), MessageEvent::Error);
	}
	Reg_Pick_Potok(passp, DBAgent::instance()->typeName(FFirstBreaks::TypeObject), m_outptFileName, "nextgenfirstbreakpicker");
	int parentId = passp->get("dataset_parent").toInt();
	QString futurePathName = DBAgentSafe::instance()->getFutureFullPath( parentId, m_outptFileName, DBAgent::instance()->typeName(FFirstBreaks::TypeObject) );
	QString xKey, yKey, azimuth;
	coordinateWords(m_key, xKey, yKey);
	outputFile.setFileName(futurePathName);
	if( !outputFile.open() ) return;
	outputFile.setFormat( FirstBreaksNameSpace::standartFormat() );

	QString futurePathNameEnergyCenters;
	QString futurePathNameUpBorder;
	QString futurePathNameDownBorder;

	if ( m_showAnalysisBorders ){
		m_outptFileNameUpBorder = m_outptFileName + "_UpBorder";
		m_outptFileNameDownBorder = m_outptFileName + "_DownBorder";
		if ( m_tuneToEnergyCenter ){
			m_outptFileNameUpBorder += "_Tuned";
			m_outptFileNameDownBorder += "_Tuned";
		}
		Reg_Pick_Potok(passp, DBAgent::instance()->typeName(FFirstBreaks::TypeObject), m_outptFileNameUpBorder, "nextgenfirstbreakpicker");
		Reg_Pick_Potok(passp, DBAgent::instance()->typeName(FFirstBreaks::TypeObject), m_outptFileNameDownBorder, "nextgenfirstbreakpicker");
		futurePathNameUpBorder = DBAgentSafe::instance()->getFutureFullPath( parentId, m_outptFileNameUpBorder, DBAgent::instance()->typeName(FFirstBreaks::TypeObject) );
		futurePathNameDownBorder = DBAgentSafe::instance()->getFutureFullPath( parentId, m_outptFileNameDownBorder, DBAgent::instance()->typeName(FFirstBreaks::TypeObject) );

		upBorderOutputFile.setFileName(futurePathNameUpBorder);
		if( !upBorderOutputFile.open() ) return;
		upBorderOutputFile.setFormat( FirstBreaksNameSpace::standartFormat() );

		downBorderOutputFile.setFileName(futurePathNameDownBorder);
		if( !downBorderOutputFile.open() ) return;
		downBorderOutputFile.setFormat( FirstBreaksNameSpace::standartFormat() );
	}
	if ( m_showEnergyCenters ){
		m_outptFileNameEnergyCenters = m_outptFileName + "_EnergyCenters";
		if ( m_tuneToEnergyCenter ){
			m_outptFileNameEnergyCenters += "_Tuned";
		}
		Reg_Pick_Potok(passp, DBAgent::instance()->typeName(FFirstBreaks::TypeObject), m_outptFileNameEnergyCenters, "nextgenfirstbreakpicker");
		futurePathNameEnergyCenters = DBAgentSafe::instance()->getFutureFullPath( parentId, m_outptFileNameEnergyCenters, DBAgent::instance()->typeName(FFirstBreaks::TypeObject) );
		energyCentersOutputFile.setFileName(futurePathNameEnergyCenters);
		if( !energyCentersOutputFile.open() ) return;
		energyCentersOutputFile.setFormat( FirstBreaksNameSpace::standartFormat() );
	}

	Trace * singleTrace;
	TraceList traceEnsemble;
	QVector<float> inputTimes;
	QMap<int, QVector<float> > outputTimes;
	outputFile.beginTransaction();
	if ( m_showAnalysisBorders ){
		upBorderOutputFile.beginTransaction();
		downBorderOutputFile.beginTransaction();
	}
	if ( m_showEnergyCenters ){
		energyCentersOutputFile.beginTransaction();
	}

	trace_buffer = new float[ m_maxTraceNum * m_nSamples ];
	WaveformAnalyzer wfAnal;
	QVector<ExtremumCurve> extremumCurveVector;
	QMap<int, QVector<FoundExtremum> > foundExtremumMap;
	QMap<int, QPointF> traceCoordinates; // [x, y, azimuth]
	QPointF coordinates;
	float *bufferPtr;
	setTraceHeaderMap( *thm );
	setPassport( *passp );
	QVector<float> buffres;

	std::vector<float>signalToAutocorr;
	std::vector<float>autocorrRes;

	while ( singleTrace= getOne() ) {

		if( isStopped() ) {
			return;
		}
		if (singleTrace->isTheLastOne()) {
			if( m_doesPut ){
				TraceHeader thLast(thm);
				thLast.headerCopy( singleTrace->getTraceHeader() );
				Trace* tout = new Trace (thLast,singleTrace->getTraceData());
				tout->setTheLastOne();
				putOne(tout);
				delOne( singleTrace );
			}
			continue;
		}
		if(singleTrace->isEnsembleEnd()) {
			int sou = static_cast<int>( singleTrace->getDouble( TraceNameSpace::toString( TraceNameSpace::SOURCE ) ) );
			coeffsMap.clear();
			traceParamsMap.clear();
			points.clear();
			groups_.clear();
			traceEnsemble.append(singleTrace);
			int traceCount = traceEnsemble.count();

			float tempTime;
			for(int i = 0; i < traceCount; ++i)
			{
				memcpy(trace_buffer + i * m_nSamples, traceEnsemble[i]->getDataPointer(), sizeof(float)*m_nSamples);
				float radius = 0.0f;
				if ( isAoffset ) radius = (int)traceEnsemble[i]->getDouble(TraceNameSpace::toString( TraceNameSpace::AOFFSET ));
				if ( !isAoffset && isOffset ) radius = abs((int)traceEnsemble[i]->getDouble(TraceNameSpace::toString( TraceNameSpace::OFFSET )));
				tempTime = rayInterpolator.interpolateLinear(traceEnsemble[i]->getDouble(m_key),
															 calcAzimuthInside(traceEnsemble[i]),
															 radius );
				inputTimes.push_back( ceil( tempTime ) );
				coordinates.ry() = radius;
				coordinates.ry() = calcAzimuthInside(traceEnsemble[i]);
				traceCoordinates.insert( i, coordinates );
				SeisTraceParams tp;
				tp.azimuth = calcAzimuthInside(traceEnsemble[i]);
				tp.dtMs = m_dt;
				tp.radius = radius;
				traceParamsMap.insert( i, tp);
			}
			QVector<float> refTimes;
			wfAnal.setWaveform( trace_buffer, traceCount, m_nSamples, m_dt, m_maxTraceNum );
			wfAnal.setParameters( m_timeWindow, m_startDelay, 0.00001f, m_maxExpectedExtremumsPerWindow, m_usePreprocessing, 3.0f, m_useEnsembleNormalization  );
			wfAnal.preprocessPhaseActivate( inputTimes );
			int lengthOfBuffer = wfAnal.getMaxLength();
			if ( buffres.capacity() < wfAnal.getMaxLength()*m_maxTraceNum ){
				buffres.reserve(wfAnal.getMaxLength()*m_maxTraceNum);
			}
			if ( signalToAutocorr.capacity() < lengthOfBuffer ){
				signalToAutocorr.reserve(lengthOfBuffer);
			}
			if ( autocorrRes.capacity() < lengthOfBuffer ){
				autocorrRes.reserve(lengthOfBuffer);
			}

			int startIndex = -1;

			QVector<QVariant> point;

			std::vector<float> startPoints;
			std::vector<float> energyCenters;
			std::vector<float> Autocorrelation_FFTDuration;

			startPoints.reserve( traceCount );
			energyCenters.reserve( traceCount );
			Autocorrelation_FFTDuration.reserve( traceCount );

			signalToAutocorr.resize(lengthOfBuffer);
			autocorrRes.resize(lengthOfBuffer);

			float rmsSignalDuration = 0;
			for (int i = 0; i < traceCount; ++i) {
				bufferPtr = wfAnal.getProcessedBufferForTrace(i, startIndex);
				doFrontSmooth(bufferPtr, lengthOfBuffer);
				doBackSmooth(bufferPtr, lengthOfBuffer);
				int tempEnergyCenter = qFloor( findApproximateEnergyCenter( bufferPtr, 0, lengthOfBuffer) ) ;
				if( startIndex<0 ){
					energyCenters.push_back( static_cast<float>(qAbs(tempEnergyCenter+startIndex)*m_dt) );
				}else{
					energyCenters.push_back( static_cast<float>((tempEnergyCenter+startIndex)*m_dt) );
				}
				memcpy( &signalToAutocorr.at(0), bufferPtr, sizeof(float) * lengthOfBuffer );
				autocorrRes = calculateAutocorrelation_FFT(signalToAutocorr);
				int tempSignalDuration = findDominantPeriod(autocorrRes,1);
				Autocorrelation_FFTDuration.push_back(static_cast<float>(tempSignalDuration));
				rmsSignalDuration += powf(static_cast<float>(tempSignalDuration*m_signalWidthCoefficent*m_dt),2);
			}
			rmsSignalDuration/=traceCount;
			rmsSignalDuration = sqrtf(rmsSignalDuration);
			send( QString( "RMS signal duration %1 at %2 : %3" ).arg(rmsSignalDuration)
																.arg( m_ensembleFirstKey )
																.arg( (int)(traceEnsemble.at(0)->getDouble( m_ensembleFirstKey )) ));
			if( m_tuneToEnergyCenter ){
				wfAnal.setWaveform( trace_buffer, traceCount, m_nSamples, m_dt, m_maxTraceNum );
				wfAnal.setParameters( rmsSignalDuration * 2 + 1, -1*rmsSignalDuration, 0.00001f, m_maxExpectedExtremumsPerWindow, m_usePreprocessing, 3.0f, m_useEnsembleNormalization );
				wfAnal.preprocessPhaseActivate( QVector<float>::fromStdVector( energyCenters ) );
				lengthOfBuffer = wfAnal.getMaxLength();
				if ( buffres.capacity() < wfAnal.getMaxLength()*m_maxTraceNum ){
					buffres.reserve(wfAnal.getMaxLength()*m_maxTraceNum);
				}
				if ( signalToAutocorr.capacity() < lengthOfBuffer ){
					signalToAutocorr.reserve(lengthOfBuffer);
				}
				if ( autocorrRes.capacity() < lengthOfBuffer ){
					autocorrRes.reserve(lengthOfBuffer);
				}
				energyCenters.clear();
			}
			for (int i = 0; i < traceCount; ++i) {
				bufferPtr = wfAnal.getProcessedBufferForTrace(i, startIndex);
				doFrontSmooth(bufferPtr, lengthOfBuffer);
				doBackSmooth(bufferPtr, lengthOfBuffer);
				int tempEnergyCenter = qFloor( findApproximateEnergyCenter( bufferPtr, 0, lengthOfBuffer) ) ;

				if( startIndex<0 ){
					energyCenters.push_back( static_cast<float>(qAbs(tempEnergyCenter+startIndex)*m_dt) );
				}else{
					energyCenters.push_back( static_cast<float>((tempEnergyCenter+startIndex)*m_dt) );
				}
				memcpy( &signalToAutocorr.at(0), bufferPtr, sizeof(float) * lengthOfBuffer );
				autocorrRes = calculateAutocorrelation_FFT(signalToAutocorr);
				int tempSignalDuration = 0;

				if ( m_tuneToEnergyCenter ){
					tempSignalDuration = qFloor( rmsSignalDuration / 2 );
				}else{
					tempSignalDuration = findDominantPeriod(autocorrRes,1);
				}
				if ( m_tuneToEnergyCenter ){
					doFrontSmooth(bufferPtr, tempEnergyCenter, tempSignalDuration);
				}else{
					doFrontSmooth(bufferPtr, tempEnergyCenter, tempSignalDuration*m_signalWidthCoefficent );
				}
				Eigen::VectorXf tempCoefs;
				if ( m_tuneToEnergyCenter )
				{
					tempCoefs = FAWSolver::fitLegendreApproximationExtended(bufferPtr, lengthOfBuffer, qFloor(tempEnergyCenter),
																			qFloor(tempSignalDuration),
																			m_polynomDegree, m_lambda, m_weightingMode,
																			m_degreePenalty, m_compressionCoefficient );
				} else{
					tempCoefs = FAWSolver::fitLegendreApproximationExtended(bufferPtr, lengthOfBuffer, qFloor(tempEnergyCenter),
																			qFloor(tempSignalDuration) * m_signalWidthCoefficent,
																			m_polynomDegree, m_lambda, m_weightingMode,
																			m_degreePenalty, m_compressionCoefficient );
				}
				coeffsMap.insert( i, tempCoefs );
				traceParamsMap[i].windowCenterTime = energyCenters.at(i);
				if( m_tuneToEnergyCenter ){
					traceParamsMap[i].windowHalfWidth = qFloor(tempSignalDuration);
				}else{
					traceParamsMap[i].windowHalfWidth = qFloor(tempSignalDuration) * m_signalWidthCoefficent;
				}
				if ( m_showApproximation ){
					memset( bufferPtr, 0, sizeof(float)*lengthOfBuffer );
					if ( m_tuneToEnergyCenter ){
						FAWSolver::synthesizeLegendreFragmentInPlace( tempCoefs, bufferPtr, lengthOfBuffer,	 qFloor(tempEnergyCenter), qFloor(tempSignalDuration) );
					} else {
						FAWSolver::synthesizeLegendreFragmentInPlace( tempCoefs, bufferPtr, lengthOfBuffer,	 qFloor(tempEnergyCenter), qFloor(tempSignalDuration)*m_signalWidthCoefficent );
					}
				}
				point.clear();
				point.push_back((int)traceEnsemble[i]->getDouble(m_key));
				point.push_back((int)traceEnsemble[i]->getDouble(xKey)/100);
				point.push_back((int)traceEnsemble[i]->getDouble(yKey)/100);
				if ( isAoffset ) point.push_back((int)traceEnsemble[i]->getDouble(TraceNameSpace::toString( TraceNameSpace::AOFFSET )));
				if ( !isAoffset && isOffset ) point.push_back(abs((int)traceEnsemble[i]->getDouble(TraceNameSpace::toString( TraceNameSpace::OFFSET ))));
				point.push_back((int)(calcAzimuthInside(traceEnsemble[i])));
				point.push_back(-1);
				//				if (outputTimes.contains(i)) point[5] = (int)(outputTimes.value(i)[0]);//*1000 ????
				point[5] = 500;
				point.push_back( (int)traceEnsemble[i]->getDouble(m_traceUniqueKeyIndex) );
				points.insert( i, point );
				//				if ( point[5].toInt()>0 )outputFile.addRelation(point);
				if ( m_showEnergyCenters ) {
					point[5] = energyCenters.at(i);
					energyCentersOutputFile.addRelation(point);
				}
				if ( m_showAnalysisBorders ) {
					if ( m_tuneToEnergyCenter ){
						point[5] = max( energyCenters.at(i) - tempSignalDuration*m_dt, 0 );
						upBorderOutputFile.addRelation(point);
						point[5] = energyCenters.at(i) + tempSignalDuration*m_dt;
						downBorderOutputFile.addRelation(point);
					}
					else{
						point[5] = max( energyCenters.at(i) - tempSignalDuration*m_signalWidthCoefficent *m_dt, 0 );
						upBorderOutputFile.addRelation(point);
						point[5] = energyCenters.at(i) + tempSignalDuration*m_signalWidthCoefficent *m_dt;
						downBorderOutputFile.addRelation(point);
					}

				}
				if( m_doesPut ){
					writeOutTrace( traceEnsemble.at(i), thm, bufferPtr, startIndex, lengthOfBuffer);
					delOne( traceEnsemble.at(i) );
				}
				//				float tTime = extremaIndex.getFirstNegExtrema(i).timeMs;
				//				point[5] = energyCenters.at(i) - 20.0f;
				//				outputFile.addRelation(points[i]);
			}
			//			extremaIndex.buildNeighborMap_( traceParamsMap, 2 );
			//			QMap<int, VelocityResult> velRes = runVelocityAnalysis( coeffsMap, traceParamsMap, extremaIndex, m_dt );

			extremaIndex.build(coeffsMap, traceParamsMap, m_maxExpectedExtremumsPerWindow, 3);
			extremaIndex.computeZeroCrossingsPrecise(coeffsMap, traceParamsMap);
			extremaIndex.computeExtremaMetricsAnalytic( coeffsMap, traceParamsMap );
			SectorSurfaceBuilder ssb(10, 25, false);
			QVector<SectorSurface> secSurf = ssb.build( coeffsMap, traceParamsMap, extremaIndex);
			groups_ = ssb.groups();
			secSurf_ = secSurf;

			foreach (QVector<int> indexesInGroup, groups_.values()) {
				QVector<int> group = indexesInGroup;

				if (group.size() <= 1) continue;

				int base = m_spatialApproximationBase;
				if (base > 1 && group.size() == base - 1) {
					base = qMax(3, (base / 2) | 1);
				}

				if (group.size() < base) {
					continue;
				}


				QVector<QPair<double,int> > rlist;
				for (int k = 0; k < indexesInGroup.size(); ++k) {
					int tid = indexesInGroup[k];
					if (!traceParamsMap.contains(tid)) continue;
					double rr = traceParamsMap.value(tid).radius;
					rlist.push_back(qMakePair(rr, tid));
				}
				if (rlist.isEmpty()) continue;

				struct CompareByRadiusFunctor {
					bool operator()(const QPair<double, int>& a, const QPair<double, int>& b) const {
						return a.first < b.first;
					}
				};

				std::sort(rlist.begin(), rlist.end(), CompareByRadiusFunctor());


				std::vector<int> traceIds;
				traceIds.reserve(rlist.size());
				for (const auto& p : rlist) {
					traceIds.push_back(p.second);
				}

				int halfWin = (base - 1) / 2;

				for (int centerIdx = 0; centerIdx < traceIds.size(); ++centerIdx) {
					int centerTraceId = traceIds[centerIdx];

					std::vector<const float*> buffers;
					std::vector<float> weights;
					buffers.reserve(base);
					weights.reserve(base);

					for (int offset = -halfWin; offset <= halfWin; ++offset) {
						int neighIdx = centerIdx + offset;
						if (neighIdx < 0) {
							int refl = centerIdx - offset;
							if (refl >= traceIds.size()) refl = 0;
              
							neighIdx = refl;
						}
						if (neighIdx >= traceIds.size()) {
							int refl = centerIdx - offset;
							if (refl < 0) refl = traceIds.size()-1;
							neighIdx = refl;
						}

						int neighTraceId = traceIds[neighIdx];
						const float* buf = wfAnal.getProcessedBufferForTrace(neighTraceId, startIndex);

						if (!buf) {
							buffers.clear();
							break;
						}
						buffers.push_back(buf);
						weights.push_back(1.0f / (1.0f + std::abs(offset)));
					}

					if (buffers.size() < base) continue;

					int bufferLength = lengthOfBuffer;
					int centerSample = qFloor(findApproximateEnergyCenter(buffers[halfWin], 0, bufferLength));

					int winHalfStacked;
					if (m_tuneToEnergyCenter) {
						winHalfStacked = qFloor(rmsSignalDuration / 2);
					} else {
						const float* centerBuffer = buffers[halfWin];
						memcpy(&signalToAutocorr.at(0), centerBuffer, sizeof(float) * bufferLength);
						autocorrRes = calculateAutocorrelation_FFT(signalToAutocorr);
						winHalfStacked = findDominantPeriod(autocorrRes, 1) * m_signalWidthCoefficent;
					}

					Eigen::VectorXf coeffs = FAWSolver::fitStackedLegendreApproximation(

								 buffers,
								 weights,
								 bufferLength,
								 centerSample,
								 winHalfStacked,
								 m_polynomDegree,
								 m_lambda,
								 m_weightingMode,
								 m_degreePenalty,
								 m_compressionCoefficient
								 );

					coeffsMap.insert(centerTraceId, coeffs);
					traceParamsMap[centerTraceId].windowHalfWidth = winHalfStacked;
				}
			}


			foreach (QVector<int> indexesInGroup, groups_.values()) {
				QVector<int> group = indexesInGroup;
				QMap<int, FoundExtremum> chosen = processAzimuthGroupOptimizeSmoothness(
													  group, traceParamsMap, extremaIndex,
													  m_extrNumber, m_pilotMode,
													  m_timeWeight, m_widthWeight,
													  m_energyWeight, m_signConsistencyWeight, m_maxIter);

				foreach (FoundExtremum fe, chosen) {
					float time = fe.timeMs;
					points[fe.traceIndex][5] = time;
					outputFile.addRelation(points[fe.traceIndex]);
				}
			}
			points.clear();
			outputTimes.clear();
			traceEnsemble.clear();
			inputTimes.clear();
			continue;
		}

		traceEnsemble.append(singleTrace);
		continue;
	} // while(singleTrace = getOne())
	outputFile.endTransaction();
	outputFile.close();
	if ( m_showAnalysisBorders ){
		upBorderOutputFile.endTransaction();
		upBorderOutputFile.close();

		downBorderOutputFile.endTransaction();
		downBorderOutputFile.close();
	}
	if ( m_showEnergyCenters  ){
		energyCentersOutputFile.endTransaction();
		energyCentersOutputFile.close();
	}
	putOne(0);
	//    delete bufferPtr;
	delete [] trace_buffer;
	send(QObject::tr("fb ensemble pick is finished"));
	GuiCore::appendLastUsedWords(m_outptFileName);
	addMsg("");

	return;
}